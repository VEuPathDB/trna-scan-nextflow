#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { tRNAScan } from './modules/tRNAScan.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------


workflow {
  // Param checking
  if(!params.fastaSubsetSize) {
    throw new Exception("Missing params.fastaSubsetSize")
  }
  if(!params.inputFilePath) {
    throw new Exception("Missing params.inputFilePath")
  }

  seqs = Channel.fromPath(params.inputFilePath)
                .splitFasta(by: params.fastaSubsetSize, file: true)

  tRNAScan(seqs)
}
