#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runtRNAScan {
  container = 'veupathdb/trnascan:1.0.0'
  
  input:
  path subsetFasta

  output:
  path (params.outputFileName), emit: tab
  path ('subset.gff'), emit: gff

  script:
    template 'runtRNAScan.bash'
}


process indexGff {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'

  publishDir params.outputDir, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gz'
    path '*.gz.tbi'

  script:
  """
  sort -k1,1 -k4,4n -k3,3r $gff > ${outputFileName}
  bgzip ${outputFileName}
  tabix -p gff ${outputFileName}.gz
  """
}


workflow tRNAScan {

  take:
  seqs

  main:
  trnascanResults = runtRNAScan(seqs)
  trnascanResults.tab.collectFile(keepHeader: true, skip: 3, storeDir: params.outputDir)
  indexGff(trnascanResults.gff.collectFile(keepHeader: false, skip: 1), params.outputGFFName);
}
