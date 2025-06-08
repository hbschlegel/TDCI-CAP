#!/usr/bin/env bash

# Exit on any errors
set -euo pipefail

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
