#!/usr/bin/env bash

# Exit on any errors
set -euo pipefail

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$(dirname "$0")/../dep/hdf5-1.14.4-3/hdf5-install/lib"
export PATH=$PATH:"$(dirname "$0")/../dep/hdf5-1.14.4-3/hdf5-install/lib"

if [ -f "data.h5" ]; then
  mv data.h5 data.1.h5
fi

# Execute TDCI
"$(dirname "$0")/tdci_core" "$@"

# Combine hdf5 files
for i in $(seq 1 $(grep -oP 'ndir\s*=\s*\K\d+' input)); do 
  h5copy -i thread$i.h5 -o data.h5 -s /direction_$i -d /direction_$i
  rm -f thread$i.h5
done

h5copy -i metadata.h5 -o data.h5 -s /metadata -d /metadata
rm -f metadata.h5

echo "TDCI Wrapper Exit"
