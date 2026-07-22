# Design: EukHighConfidenceFilter-based tRNA filtering

Date: 2026-07-22
Branch: `test_trnascan_filter`
Status: approved for planning

## Problem

The current `strictFilter` process reimplements tRNAscan-SE's high-confidence
filtering with a hand-rolled `awk` expression. It filters the Note column for
strings (`isotype mismatch`, `unexpected anticodon`) that are written **only by
`EukHighConfidenceFilter`** — never by a plain `tRNAscan-SE` run. On the raw
output this pipeline produces, those two criteria are inert, so the filter is
weaker than it appears and does not match the high-confidence definition used
for genomes submitted to NCBI / published in GtRNAdb.

At the same time the pipeline runs `tRNAscan-SE` with `-H --detail -f`, which
computes exactly the isotype-specific and secondary-structure scores that
`EukHighConfidenceFilter` needs — and then discards them (the `.ss` file is
emitted and never consumed).

## Goal

Replace the ad-hoc strict filter with the bundled `EukHighConfidenceFilter`
(shipped in the `veupathdb/trnascan:1.0.0` container at
`/usr/bin/EukHighConfidenceFilter`), retain only its `high confidence set`
category, and validate the retained counts against local gold-standard genomes.

## Decisions (locked during brainstorming)

- **Retention policy:** publish only tRNAs tagged `high confidence set`.
- **Thresholds:** exposed as nextflow params, defaulting to the tool defaults
  (`cmscore1=50`, `ssscore1=10`, `isoscore1=70`).
- **Invocation mode:** tag-only (no `-r`), then select the passing set
  ourselves, so the full tagged output and the `.log` category breakdown remain
  available for diagnosing count gaps.
- **Comparison metric:** total retained tRNA count per genome vs the gold count.
- **Toggle:** keep `simpleFilter` as the `applyHighConfFilter=false` path; flip
  the default to `true`.

## EukHighConfidenceFilter interface (verified from container)

```
EukHighConfidenceFilter -i <result.out> -s <ss file> -o <out dir> -p <prefix>
  -c1 <domain/overall score cutoff>   default 50
  -m1 <secondary structure cutoff>    default 10
  -e1 <isotype-specific cutoff>       default 70
  -r                                  remove instead of tag (NOT used)
```
Requires the tabular `.out` and the `.ss` file **together** (cross-keyed by
sequence name), so it must run on the merged result, not per chunk. Produces
`<prefix>.out`, `<prefix>.ss`, `<prefix>.log`.

## Process changes

| Process | Change |
|---|---|
| `hardMask`, `runtRNAScan`, `mergeTab`, `indexGff` | unchanged |
| `mergeSs` | **new** — concatenate per-chunk `subset.ss` into `merged.ss` |
| `mergeGff` | **new process** replacing the `collectFile(name:'merged.gff', skip:1)` call; strip the full `##` header block, emit exactly one `##gff-version 3` |
| `eukHighConfidenceFilter` | **new** — replaces `strictFilter` |
| `simpleFilter` | kept for the `applyHighConfFilter=false` path |
| `filterGff` | **new shared process** — coord-lookup + GFF join, extracted so strict and simple paths call it once instead of duplicating the awk |

### `eukHighConfidenceFilter`

```bash
EukHighConfidenceFilter -i ${mergedTab} -s ${mergedSs} -o . -p hiConf \
    -c1 ${params.cmScore} -m1 ${params.ssScore} -e1 ${params.isoScore}

# Select the passing set. Match is LINE-LEVEL, not $NF: the tag
# "high confidence set" is three words, so $NF would only see "set".
head -3 hiConf.out > tRNAScan.out
awk 'NR>3 && /high confidence set/' hiConf.out >> tRNAScan.out

echo "Retained $(awk 'NR>3' tRNAScan.out | wc -l) tRNAs (high confidence set)" >&2
```

Outputs: `tRNAScan.out` (published), `hiConf.log` (published — per-category
counts for diagnosing gaps), and the coordinate list handed to `filterGff`.

### `filterGff` (shared)

Takes the retained `.out` and the merged GFF; builds the strand-normalized
`seq/min/max` coordinate lookup and emits `filtered.gff`. Both `simpleFilter`
and `eukHighConfidenceFilter` feed it, removing the duplicated join logic.

## New params (`nextflow.config`)

```
applyHighConfFilter = true      // flipped from false
cmScore  = 50                   // EukHCF -c1 (overall/domain score)
ssScore  = 10                   // EukHCF -m1 (secondary structure)
isoScore = 70                   // EukHCF -e1 (isotype-specific)
```
`minInfScore` is retained but now used only by the `simpleFilter` path.

## Data flow

```
hardMask
  -> runtRNAScan
       -> { mergeTab, mergeSs, mergeGff }
            -> eukHighConfidenceFilter (applyHighConfFilter=true)
               |  or simpleFilter       (applyHighConfFilter=false)
               -> filterGff
                    -> indexGff
```

## Comparison harness (validation, not committed pipeline code)

A standalone script kept in the scratchpad:

1. Input: a two-column list of `gold_genome_fasta  expected_count` (paths on
   this machine, provided by the user).
2. For each genome: run the pipeline with `applyHighConfFilter=true`, count
   `awk 'NR>3' tRNAScan.out | wc -l`.
3. Emit a table: `genome | pipeline_count | gold_count | delta`, plus the
   EukHCF `.log` category counts (how many dropped as pseudo /
   secondary-filtered / mismatch / unexpected-anticodon).

Blocking input for the run (not for the spec): gold FASTA paths + expected
counts.

## Documentation / cleanup

- Update the stale `README.md` mermaid diagram and param table to reflect the
  hardMask → merge → filter → index flow and the new params.
- Delete `modules/templates/fixHeader.bash` (no longer referenced).

## Risks / to verify during implementation

- `mergeTab` reproduces the tRNAscan-SE header by hand; confirm
  `EukHighConfidenceFilter` parses that merged `.out` header correctly.
- Confirm concatenated `merged.ss` records stay correctly keyed to the merged
  `.out` (sequence-name + coordinate match) across chunk boundaries.
- Confirm the retained-count of `high confidence set` on a known genome is in a
  sane range before trusting the gold comparison.

## Out of scope

- Coordinate-level (TP/FP/FN) concordance against gold sets — counts only for
  now; concordance is a possible follow-up.
- Non-eukaryotic / bacterial/archaeal filtering modes.
