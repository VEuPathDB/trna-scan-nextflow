nextflow.enable.dsl=2

process tRNAscan {
    input:
    path 'subset.fa'
    output:
    path 'subset.scanned'
    """
    tRNAscan-SE subset.fa -o subset.scanned 
    sed -i '1,3d' subset.scanned
    """
}

process addHeader {
  publishDir params.outputDir, mode: 'copy', saveAs: {filename -> params.outputFileName}
  input:
  path 'subset.scanned_noheader'
  output:
  path 'header.txt'
  """
  echo -e "Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron\tBounds\tInf\nName\t\ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\t\n--------\t------\t-----\t------\t----\t-----\t-----\t----\t------\t------" | cat - subset.scanned_noheader > header.txt
  """
}

workflow {
  seqs = channel.fromPath(params.inputFilePath).splitFasta(by:params.numberOfSplits, file:true)
  tRNAscan(seqs).collectFile(name: params.outputFileName) | addHeader
}

