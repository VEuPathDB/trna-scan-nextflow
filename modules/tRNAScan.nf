#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runtRNAScan {
  container = 'veupathdb/trnascan:latest'
  
  input:
    path subsetFasta

  output:
    path 'subset.scanned'

  script:
    template 'runtRNAScan.bash'
}


process fixHeader {
  container = 'veupathdb/trnascan:latest'

  publishDir params.outputDir, saveAs: {filename->params.outputFile}, mode: 'copy'

  input:
    path fileWithHeader

  output:
    path 'datafile'

  script:
    template 'fixHeader.bash'
}


workflow tRNAScan {

  take:
    seqs

  main:

    runtRNAScan(seqs) \
        | collectFile() \
        | fixHeader
}
