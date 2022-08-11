nextflow.enable.dsl=2

process tRNAscan {
  input:
    path 'subset.fa'

  output:
    path 'subset.scanned'

  script:
    """
    tRNAscan-SE subset.fa \
      -o subset.scanned 
    """
}

process fixHeader {
  publishDir params.outputDir, saveAs: {filename->params.outputFile}

  input:
    path 'fileWithHeader.txt'

  output:
    path 'datafile'

  script:
    """
    fixheader.pl
    """
}


workflow {
  channel.fromPath(params.inputFilePath)
    .splitFasta(by:params.fastaSubsetSize, file:true) \
      | tRNAscan \
      | collectFile() \
      | fixHeader
}

