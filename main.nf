nextflow.enable.dsl=2

process tRNAscan {
    publishDir params.outputDir, mode: 'copy', saveAs: { filename -> params.outputFileName }
    input:
    path 'subset.fa'
    output:
    path 'subset.scanned'
    """
    tRNAscan-SE subset.fa -o subset.scanned 
    """
}

workflow {
  channel.fromPath(params.inputFilePath).splitFasta(by:1, file:true) | tRNAscan
}
