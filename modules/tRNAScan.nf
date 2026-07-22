#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
// Convert soft-masked sequence to hard-masked sequence
// Input genomes are already soft-masked (repeats in lowercase); tRNAscan-SE
// ignores case, so soft masking alone has no effect. Converting lowercase to N
// hard-masks those regions, keeping tRNA-derived SINE tandem arrays from
// flooding the output on large repeat-rich genomes.
//
// This is a pure per-record transform, so it runs on each FASTA chunk
// independently (no need to see the whole genome at once).
// ---------------------------------------------------------------------------
process hardMask {
  container 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

  input:
  path subsetFasta

  output:
  path 'hardmasked.fa', emit: masked

  script:
  """
  seqkit replace -s -p '[a-z]' -r N -w 0 ${subsetFasta} > hardmasked.fa
  """
}

// ---------------------------------------------------------------------------
// Run tRNAscan-SE on each FASTA chunk
// -E  : eukaryotic mode
// -H  : include HMM/ss breakdown scores (required by EukHighConfidenceFilter
//       and by our strict post-filter which checks col12 isotype)
// --detail : include isotype-specific model scores (required for col12/col13)
// -f  : secondary structure output (required by EukHighConfidenceFilter)
// ---------------------------------------------------------------------------
process runtRNAScan {
  container 'veupathdb/trnascan:1.0.0'

  input:
  path subsetFasta

  output:
  path (params.outputFileName), emit: tab
  path ('subset.gff'),          emit: gff
  path ('subset.ss'),           emit: ss

  script:
    template 'runtRNAScan.bash'
}

// ---------------------------------------------------------------------------
// Merge per-chunk tabular outputs with an explicit known-good header
// Avoids the collectFile(keepHeader:true) bug where an empty-hit chunk
// can corrupt or lose the header (seen with A. gambiae PEST)
// ---------------------------------------------------------------------------
process mergeTab {
  input:
  path tabFiles, stageAs: '*.tab'

  output:
  path 'merged.tab'

  script:
  """
  # Take the 3-line header from the first chunk so the column layout matches
  # tRNAscan-SE's real output exactly. A hand-written header previously dropped
  # the empty field in row 1 (real row1 has 13 tab-separated fields, an empty
  # field 2 from the double-tab after "Sequence"), which shifted the Anti/score
  # columns and broke EukHighConfidenceFilter's column parsing.
  # Every chunk emits the full 3-line header even with zero hits, so the first
  # chunk is a safe header source; only data rows (line 4+) are appended.
  first=1
  for f in ${tabFiles}; do
    if [ "\$first" -eq 1 ]; then
      head -3 \$f > merged.tab
      first=0
    fi
    tail -n +4 \$f >> merged.tab
  done
  """
}

// ---------------------------------------------------------------------------
// Concatenate per-chunk secondary-structure files into one .ss file.
// EukHighConfidenceFilter cross-keys .ss records (by sequence name) against the
// merged tabular output, so it needs the whole set in one file.
// ---------------------------------------------------------------------------
process mergeSs {
  input:
  path ssFiles, stageAs: '*.ss'

  output:
  path 'merged.ss'

  script:
  """
  cat ${ssFiles} > merged.ss
  """
}

// ---------------------------------------------------------------------------
// Merge per-chunk GFF3 with a single clean header.
// Replaces collectFile(skip:1), which stripped only one line per chunk and left
// duplicate/renumbered ## header lines interleaved in the output.
// ---------------------------------------------------------------------------
process mergeGff {
  input:
  path gffFiles, stageAs: '*.gff'

  output:
  path 'merged.gff'

  script:
  """
  echo '##gff-version 3' > merged.gff
  for f in ${gffFiles}; do
    grep -v '^#' \$f >> merged.gff || true
  done
  """
}

// ---------------------------------------------------------------------------
// High-confidence filter using tRNAscan-SE's bundled EukHighConfidenceFilter.
// Run in tag-only mode (no -r) so the full tagged output and the .log category
// breakdown survive for diagnosing count gaps; then retain only the
// "high confidence set". The tag is three words, so the match is LINE-LEVEL,
// not $NF (which would only see "set").
// ---------------------------------------------------------------------------
process eukHighConfidenceFilter {
  container 'veupathdb/trnascan:1.0.0'
  publishDir params.outputDir, mode: 'copy', pattern: 'tRNAScan.out'
  publishDir params.outputDir, mode: 'copy', pattern: 'hiConf.log'

  input:
  path mergedTab
  path mergedSs

  output:
  path 'tRNAScan.out', emit: tab
  path 'hiConf.log',   emit: log

  script:
  """
  EukHighConfidenceFilter -i ${mergedTab} -s ${mergedSs} -o . -p hiConf \\
      -c1 ${params.cmScore} -m1 ${params.ssScore} -e1 ${params.isoScore}

  head -3 hiConf.out > tRNAScan.out
  awk 'NR>3 && /high confidence set/' hiConf.out >> tRNAScan.out

  echo "Retained \$(awk 'NR>3' tRNAScan.out | wc -l) tRNAs (high confidence set)" >&2
  """
}

// ---------------------------------------------------------------------------
// Apply a strict flag-based filter directly on the tabular output:
//
//   col5  = Type (anticodon-based isotype call)
//   col9  = Inf Score
//   col12 = Isotype (model-based isotype call)  [shifted by -H vs old format]
//   colNF = Note (last field - variable position with --detail -H)
//
// Kept entries must satisfy ALL of:
//   1. Not pseudo          (colNF !~ /pseudo/)
//   2. No isotype mismatch (colNF !~ /isotype.mismatch/)
//   3. No unexpected anticodon (colNF !~ /unexpected.anticodon/)
//   4. No truncation       (colNF !~ /trunc/)
//   5. Inf score >= params.minInfScore (default 60)
//   6. Anticodon-based type matches isotype model (col5 == col12)
//
// This combination ensures only genuinely high-quality, unambiguous tRNA
// gene predictions are retained, regardless of genome size or repeat content.
// After RepeatMasking, real tRNAs will pass all these criteria cleanly.
// ---------------------------------------------------------------------------
process strictFilter {
  container 'veupathdb/trnascan:1.0.0'
  publishDir params.outputDir, mode: 'copy', pattern: 'tRNAScan.out'

  input:
  path mergedTab
  path mergedGff

  output:
  path 'tRNAScan.out', emit: tab
  path 'filtered.gff', emit: gff

  script:
  """
  # Write header
  head -3 ${mergedTab} > tRNAScan.out

  # Apply strict filters:
  # - no pseudo, no isotype mismatch, no unexpected anticodon, no truncation
  # - Inf score >= minInfScore
  # - anticodon-based type (col5) must match isotype model (col12)
  awk -v minscore=${params.minInfScore} '
    NR>3 &&
    \$NF !~ /pseudo/ &&
    \$NF !~ /isotype.mismatch/ &&
    \$NF !~ /unexpected.anticodon/ &&
    \$NF !~ /trunc/ &&
    \$9+0 >= minscore &&
    \$5 == \$12 {print}
  ' ${mergedTab} >> tRNAScan.out

  echo "Retained \$(awk 'NR>3' tRNAScan.out | wc -l) tRNAs after strict filtering" >&2

  # Build coordinate lookup (normalize strand: always min..max to match GFF)
  awk 'NR>3 {
    start = (\$3 < \$4) ? \$3 : \$4
    end   = (\$3 < \$4) ? \$4 : \$3
    print \$1"\t"start"\t"end
  }' tRNAScan.out > hc_coords.txt

  # Filter GFF to matching coordinates
  awk 'NR==FNR {coords[\$1"\t"\$2"\t"\$3]=1; next}
       /^#/ {print; next}
       {if (coords[\$1"\t"\$4"\t"\$5]) print}' hc_coords.txt ${mergedGff} > filtered.gff
  """
}

// ---------------------------------------------------------------------------
// Simple filter for small genomes (applyHighConfFilter = false):
// drop pseudogenes and low-scoring hits. GFF join is delegated to filterGff.
// Uses $NF for the note field to handle the --detail -H column layout.
// ---------------------------------------------------------------------------
process simpleFilter {
  container 'veupathdb/trnascan:1.0.0'
  publishDir params.outputDir, mode: 'copy', pattern: 'tRNAScan.out'

  input:
  path mergedTab

  output:
  path 'tRNAScan.out', emit: tab

  script:
  """
  head -3 ${mergedTab} > tRNAScan.out
  awk -v minscore=${params.minInfScore} '
    NR>3 &&
    \$NF !~ /pseudo/ &&
    \$9+0 >= minscore {print}
  ' ${mergedTab} >> tRNAScan.out

  echo "Retained \$(awk 'NR>3' tRNAScan.out | wc -l) tRNAs after simple filtering" >&2
  """
}

// ---------------------------------------------------------------------------
// Filter the merged GFF to the coordinates retained in the filtered .out.
// Shared by both the strict (EukHCF) and simple filter paths.
// Strand is normalized to min..max to match GFF start<end convention.
// ---------------------------------------------------------------------------
process filterGff {
  publishDir params.outputDir, mode: 'copy', pattern: 'filtered.gff'

  input:
  path retainedTab
  path mergedGff

  output:
  path 'filtered.gff'

  script:
  """
  awk 'NR==FNR { if (FNR>3) ids[\$1".trna"\$2]=1; next }
       /^#/ { print; next }
       { if (match(\$9, /[A-Za-z0-9_.]+[.]trna[0-9]+/) && (substr(\$9,RSTART,RLENGTH) in ids)) print }' ${retainedTab} ${mergedGff} > filtered.gff
  """
}

// ---------------------------------------------------------------------------
// Sort, bgzip and tabix-index the filtered GFF
// ---------------------------------------------------------------------------
process indexGff {
  container 'biocontainers/tabix:v1.9-11-deb_cv1'
  publishDir params.outputDir, mode: 'copy'

  input:
  path gff
  val outputFileName

  output:
  path '*.gz'
  path '*.gz.tbi'

  script:
  """
  sort -k1,1 -k4,4n -k3,3r ${gff} > ${outputFileName}
  bgzip ${outputFileName}
  tabix -p gff ${outputFileName}.gz
  """
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow tRNAScan {
  take:
  seqs

  main:

  // Step 1: Convert soft mask to hard mask on each chunk
  // Repeats arrive lowercase (soft-masked); tRNAscan-SE ignores case, so we
  // convert lowercase to N to actually hard-mask them before scanning.
  maskedSeqs = hardMask(seqs)

  // Step 2: Run tRNAscan-SE on each chunk
  trnascanResults = runtRNAScan(maskedSeqs.masked)

  // Step 3: Merge outputs (tab + ss + gff), each with correct header handling
  mergedTab = mergeTab(trnascanResults.tab.collect())
  mergedSs  = mergeSs(trnascanResults.ss.collect())
  mergedGff = mergeGff(trnascanResults.gff.collect())

  // Step 4: Filter (strict = EukHighConfidenceFilter, simple = score/pseudo)
  if (params.applyHighConfFilter) {
    filtered = eukHighConfidenceFilter(mergedTab, mergedSs)
  } else {
    filtered = simpleFilter(mergedTab)
  }

  // Step 5: Sync GFF to the retained set, then index
  filteredGff = filterGff(filtered.tab, mergedGff)
  indexGff(filteredGff, params.outputGFFName)
}
