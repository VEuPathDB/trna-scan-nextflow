nextflow.enable.dsl=1
seq_qch = Channel.fromPath(params.inputFilePath).splitFasta( by:1, file:true  )

process tRNAscan {
    
    input:
    path 'subset.fa' from seq_qch
    output:
    path 'subset.fa.scanned' into scanned_qch
        
    """
    tRNAscan-SE subset.fa 
    """
}

results = scanned_qch.collectFile(storeDir: params.outputDir, name: params.outputFileName)