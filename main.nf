nextflow.enable.dsl=2

process tRNAscan {
    input:
    path 'subset.fa'
    output:
    path 'subset.scanned'
    """
    tRNAscan-SE subset.fa -o subset.scanned 
    """
}

workflow {
  channel.fromPath(params.inputFilePath).splitFasta(by:params.fastaSubsetSize, file:true) | tRNAscan | collectFile(name: params.outputFileName, skip: 3, keepHeader: true, storeDir: params.outputDir)
}

