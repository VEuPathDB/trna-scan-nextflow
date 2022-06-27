nextflow.enable.dsl=2

process tRNAscan {
    publishDir params.outputDir, mode: 'copy'
    input:
    tuple val(id), val (sequence) 
    output:
    path '*.txt'
    """
    echo -e '>$id\n$sequence' > subset.fa
    tRNAscan-SE subset.fa -o "$id"_scanned.txt 
    """
}

workflow {
  channel.fromPath(params.inputFilePath).splitFasta(by:1, record: [id: true, sequence: true]) | tRNAscan
}