# trna-scan-nextflow

A Nextflow pipeline that scans genomic sequences for transfer RNA (tRNA) genes using tRNAscan-SE.

## Overview

This pipeline identifies tRNA genes in a genome assembly by running [tRNAscan-SE 2.0](http://trna.ucsc.edu/) over the input sequences and producing a consolidated results table, a filtered GFF3 of retained tRNA loci, and an indexed, bgzipped version of that GFF3. It is used within VEuPathDB's genome annotation workflows to generate tRNA gene predictions as part of building out a genome's structural annotation. The input FASTA is split into subsets so that scanning can be parallelized across many concurrent jobs, and the individual results are merged and then filtered down to a confident tRNA set.

The pipeline supports two independent, optional steps around the core scan:

- **Hard-masking** — soft-masked (lowercase) bases can optionally be converted to `N` before scanning, to exclude those regions from tRNAscan-SE entirely. Off by default, since tRNAscan-SE ignores case and the high-confidence filter is the primary defense against tRNA-derived repeat noise.
- **Filtering** — after scanning, results are reduced to a confident tRNA set using one of two methods: tRNAscan-SE's bundled `EukHighConfidenceFilter` (strict, score/structure-based; the default), or a simpler filter that drops pseudogenes and low-scoring hits by Infernal score cutoff (useful for small genomes).

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- A container engine: Docker or Singularity (profiles are provided for both)
- Optionally, an LSF cluster (an `lsf` config profile is included for job submission)

The pipeline uses the `veupathdb/trnascan:1.0.0` container image for tRNAscan-SE and `EukHighConfidenceFilter`, a `quay.io/biocontainers/seqkit` image for hard-masking, and a `biocontainers/tabix` image for indexing the final GFF3.

## Usage

```
nextflow run VEuPathDB/trna-scan-nextflow \
  -r main \
  -entry tRNAScan \
  --inputFilePath /path/to/genomicSeqs.fa \
  --outputDir /path/to/output \
  --fastaSubsetSize 25 \
  -profile docker \
  -resume
```

The pipeline has a single workflow entry point:

- **`tRNAScan`** (default) — splits the input FASTA into subsets, optionally hard-masks each subset, runs `tRNAscan-SE` on each subset, merges the per-subset tabular/secondary-structure/GFF results, filters down to a confident tRNA set, and produces a sorted, bgzipped, tabix-indexed GFF3 of the retained loci.

## Key Parameters

| Parameter | Description |
| --- | --- |
| `params.inputFilePath` | Path to the input FASTA file of genomic sequences to scan. |
| `params.outputDir` | Directory where the final output files are published. |
| `params.outputFileName` | File name for the published tabular tRNA report (default: `output_scanned.txt`). |
| `params.outputGFFName` | File name for the published, indexed GFF3 of retained tRNA loci (default: `output_scanned.gff`). |
| `params.fastaSubsetSize` | Number of sequences to include in each FASTA subset sent to a single `tRNAscan-SE` process; controls the degree of parallelism. |
| `params.applyHardMask` | If `true`, converts soft-masked (lowercase) bases to `N` before scanning. Default `false`. |
| `params.applyHighConfFilter` | If `true` (default), filters with `EukHighConfidenceFilter` (retains only the "high confidence set"). If `false`, uses the simpler score/pseudogene filter instead — better suited to small genomes. |
| `params.minInfScore` | Infernal score cutoff used by the simple filter path when `applyHighConfFilter` is `false` (default `60`). |
| `params.cmScore` | `EukHighConfidenceFilter` domain/overall model score cutoff, `-c1` (default `50`). |
| `params.ssScore` | `EukHighConfidenceFilter` secondary structure score cutoff, `-m1` (default `10`). |
| `params.isoScore` | `EukHighConfidenceFilter` isotype-specific model score cutoff, `-e1` (default `70`). |

## Output

Published to `params.outputDir`:

- **`tRNAScan.out`** (name set by `params.outputFileName`) — a tab-delimited report of every retained tRNA, with columns for sequence name, tRNA number, genomic bounds, tRNA type, anticodon, intron bounds, and the Infernal (Inf) confidence score.
- **`hiConf.log`** — emitted only when `applyHighConfFilter` is `true`; the `EukHighConfidenceFilter` category breakdown log, useful for diagnosing filtered-out counts.
- **`<outputGFFName>.gz`** and **`<outputGFFName>.gz.tbi`** — a sorted, bgzipped GFF3 of the retained tRNA loci and its tabix index.
