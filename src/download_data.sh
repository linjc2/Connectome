#!/usr/bin/env bash
# Download the raw source data into data/raw/ and check each file's SHA-256.
# Usage (from the repository root): bash src/download_data.sh
set -euo pipefail

OPENWORM=https://raw.githubusercontent.com/openworm/ConnectomeToolbox/ffe570268b7cffb9329cd009a1bb37d60db2d465/cect/data
VARSHNEY=https://raw.githubusercontent.com/lrvarshney/elegans/07d10c2d43b3da9a1fefed2e8658a98dc4310765

fetch() {  # fetch <url> <destination> <sha256>
    if [ ! -f "$2" ]; then
        echo "Downloading $2"
        mkdir -p "$(dirname "$2")"
        curl -fsSL -o "$2" "$1"
    fi
    echo "$3  $2" | shasum -a 256 -c -
}

fetch "$OPENWORM/SI%205%20Connectome%20adjacency%20matrices.xlsx" \
      data/raw/connectome/cook2019_SI5_adjacency_matrices_original.xlsx \
      559989daa02cd9a76e9266537a6b80bfd47d338cb5d397d11296288629df364a
fetch "$OPENWORM/SI%205%20Connectome%20adjacency%20matrices%2C%20corrected%20July%202020.xlsx" \
      data/raw/connectome/cook2019_SI5_adjacency_matrices_corrected_2020-07.xlsx \
      1f4fdbf84746b69b49a8da0816f52787860ce349b638dce37924ba80f90c70c9
fetch "$OPENWORM/herm_full_edgelist.csv" \
      data/raw/connectome/herm_full_edgelist.csv \
      142693f17556148d7f962835b18ac6dd5af18b7467eef61815ebc1dd5474c0ca
fetch "$VARSHNEY/ConnOrdered_040903.mat" \
      data/raw/connectome/ConnOrdered_040903.mat \
      bd05af64c7d8c14ff44988e0f009b45d1ebe7db0b7194fff341fc4885cc91d38
fetch "$VARSHNEY/NeuronTypeOrdered_040903.mat" \
      data/raw/connectome/NeuronTypeOrdered_040903.mat \
      5baec9fa9ffad32eec0429c962f5d3e1ebb2da35ce612f91499da1dad5c03c17
fetch "https://www.wormatlas.org/papers/Altun%202009.pdf" \
      data/raw/expression/Altun2009_DevDyn_innexin_map.pdf \
      89992228914ad0c429bcea597456009d5199df9ccd25bc20c85d2ea1a30981ab
