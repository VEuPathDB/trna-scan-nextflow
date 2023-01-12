#!/usr/bin/env bash

set -euo pipefail
tRNAscan-SE $subsetFasta \
  -o subset.scanned 
