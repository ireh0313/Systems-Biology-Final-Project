#!/usr/bin/env bash
set -euo pipefail

sudo apt update
sudo apt install -y git build-essential cmake libxml2-dev libboost-all-dev libopenmpi-dev openmpi-bin

if [ ! -d "$HOME/neoParSA-master" ]; then
    git clone https://github.com/kennethabarr/neoParSA "$HOME/neoParSA-master"
fi

if [ ! -d "$HOME/transcpp-master" ]; then
    git clone https://github.com/kennethabarr/transcpp "$HOME/transcpp-master"
fi

cd "$HOME/neoParSA-master"
rm -rf build
mkdir -p build
cd build
cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_POLICY_VERSION_MINIMUM=3.5 ..
make parsa -j"$(nproc)"

cd "$HOME/transcpp-master"
make clean || true
make transcpp unfold scramble -j"$(nproc)"

ls -lh "$HOME/transcpp-master/transcpp" "$HOME/transcpp-master/unfold" "$HOME/transcpp-master/scramble"
