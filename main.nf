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

  tRNAScan(seqs)

}
