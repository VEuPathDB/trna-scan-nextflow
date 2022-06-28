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
  seqs = channel.fromPath(params.inputFilePath).splitFasta(by:params.numberOfSplits, file:true)
  tRNAscan(seqs).collectFile(name: params.outputFileName, skip: 3, keepHeader: true, storeDir: params.outputDir)
}

