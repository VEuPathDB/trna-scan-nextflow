#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
// Optional RepeatMasker hard-masking step
// Only runs when params.applyRepeatMask = true
// Uses dfam/tetools container which bundles RepeatMasker + full Dfam library
// Hard-masking (N) eliminates tRNA-derived SINE tandem arrays before tRNAscan-SE
// sees them, which is the only reliable way to handle large repeat-rich genomes
// like A. americanum (~70-80% repeat content)
// ---------------------------------------------------------------------------
process repeatMask {
  container 'dfam/tetools:latest'

  input:
  path genome

  output:
  path "${genome}.masked", emit: masked

  script:
  """
  RepeatMasker \
    -species "${params.repeatMaskSpecies}" \
    -pa 8 \
    -nolow \
    -dir . \
    ${genome}

  # If no repeats were found RepeatMasker does not produce a .masked file
  # Fall back to original so downstream processes always have input
  if [ ! -f ${genome}.masked ]; then
    cp ${genome} ${genome}.masked
  fi
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
  cat <<'EOF' > merged.tab
Sequence\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\tHMM\t2'Str\tIsotype\tIsotype\t
Name\ttRNA#\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tScore\tScore\tCM\tScore\tNote
--------\t------\t-----\t------\t----\t-----\t-----\t----\t------\t-----\t-----\t-------\t-------\t------
EOF

  for f in ${tabFiles}; do
    tail -n +4 \$f >> merged.tab
  done
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
  publishDir params.outputDir, mode: 'copy'

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
// Simple filter for small genomes (applyHighConfFilter = false)
// Still uses $NF for the note field to handle --detail -H column layout.
// ---------------------------------------------------------------------------
process simpleFilter {
  container 'veupathdb/trnascan:1.0.0'
  publishDir params.outputDir, mode: 'copy'

  input:
  path mergedTab
  path mergedGff

  output:
  path 'tRNAScan.out', emit: tab
  path 'filtered.gff', emit: gff

  script:
  """
  head -3 ${mergedTab} > tRNAScan.out
  awk -v minscore=${params.minInfScore} '
    NR>3 &&
    \$NF !~ /pseudo/ &&
    \$9+0 >= minscore {print}
  ' ${mergedTab} >> tRNAScan.out

  echo "Retained \$(awk 'NR>3' tRNAScan.out | wc -l) tRNAs after simple filtering" >&2

  awk 'NR>3 {
    start = (\$3 < \$4) ? \$3 : \$4
    end   = (\$3 < \$4) ? \$4 : \$3
    print \$1"\t"start"\t"end
  }' tRNAScan.out > hc_coords.txt

  awk 'NR==FNR {coords[\$1"\t"\$2"\t"\$3]=1; next}
       /^#/ {print; next}
       {if (coords[\$1"\t"\$4"\t"\$5]) print}' hc_coords.txt ${mergedGff} > filtered.gff
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

  // Step 1: Optional RepeatMasking
  // For large repeat-rich genomes (e.g. A. americanum, 3 Gb, ~75% repeats),
  // hard-masking SINEs before tRNAscan-SE is the only reliable way to prevent
  // tRNA-derived SINE tandem arrays from flooding the output.
  // For small genomes this step is skipped entirely.
  if (params.applyRepeatMask) {
    // seqs is the whole genome file (unsplit) - mask it first, then split
    maskedGenome = repeatMask(seqs)
    chunkedSeqs  = maskedGenome.masked
                    .splitFasta(by: params.fastaSubsetSize, file: true)
  } else {
    chunkedSeqs = seqs
  }

  // Step 2: Run tRNAscan-SE on each chunk
  trnascanResults = runtRNAScan(chunkedSeqs)

  // Step 3: Merge outputs
  // mergeTab uses explicit header to avoid collectFile(keepHeader:true) bug
  mergedTab = mergeTab(trnascanResults.tab.collect())
  mergedSs  = trnascanResults.ss.collectFile(name: 'merged.ss',  keepHeader: false, skip: 0)
  mergedGff = trnascanResults.gff.collectFile(name: 'merged.gff', keepHeader: false, skip: 1)

  // Step 4: Filter
  if (params.applyHighConfFilter) {
    filteredResults = strictFilter(mergedTab, mergedGff)
  } else {
    filteredResults = simpleFilter(mergedTab, mergedGff)
  }

  // Step 5: Index the filtered GFF (consistent with tRNAScan.out)
  indexGff(filteredResults.gff, params.outputGFFName)
}
