#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.inputFilePath) {
  seqs = Channel.fromPath(params.inputFilePath)
    .splitFasta(by:params.fastaSubsetSize, file:true)
}
else {
  throw new Exception("Missing params.inputFilePath")
}


//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { tRNAScan } from './modules/tRNAScan.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------


workflow {
  if (params.applyRepeatMask) {
    // Pass the whole genome as a single file channel - splitting happens AFTER masking
    seqs = Channel.fromPath(params.inputFilePath)
  } else {
    // Pre-split as before
    seqs = Channel.fromPath(params.inputFilePath)
                  .splitFasta(by: params.fastaSubsetSize, file: true)
  }

  tRNAScan(seqs)

}
