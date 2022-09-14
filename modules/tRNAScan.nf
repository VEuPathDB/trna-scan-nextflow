#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process runtRNAScan {
  input:
    path 'subset.fa'

  output:
    path 'subset.scanned'

  script:
    template 'runtRNAScan.bash'
}


process fixHeader {
  publishDir params.outputDir, saveAs: {filename->params.outputFile}

  input:
    path 'fileWithHeader.txt'

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