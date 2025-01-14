#!/usr/bin/env bash

set -euo pipefail

mkdir ./trnascanTmp

if [ -f "/usr/bin/tRNAscan-SE.conf" ]; then
    perl -e 'while(<>) {if(/^temp_dir:/) {print "temp_dir: trnascanTmp\n"} else {print}}' /usr/bin/tRNAscan-SE.conf  >tRNAscan-SE_mod.conf
    tRNAscan-SE $subsetFasta -o subset.scanned --conf tRNAscan-SE_mod.conf
else
  echo "ERROR:  /usr/bin/tRNAScan-SE.conf file does not exist"
  exit 1
fi
