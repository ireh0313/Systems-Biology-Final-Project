# 01. Build neoParSA and transcpp

This step rebuilds the old GitHub code on a current Ubuntu/WSL system.

Original repositories:

- `https://github.com/kennethabarr/neoParSA`
- `https://github.com/kennethabarr/transcpp`

Main compatibility fixes used during this project:

- New CMake versions reject the old `cmake_minimum_required` policy unless `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is provided.
- `neoParSA` requires MPI headers/libraries, so `libopenmpi-dev` is installed.
- `transcpp` needs `libxml2` and points its Makefile to the local `neoParSA` directory.

Run:

```bash
bash build_neoparsa_transcpp_ubuntu.sh
```

Successful build produces:

```text
~/neoParSA-master/build/lib/libparsa.a
~/transcpp-master/transcpp
~/transcpp-master/unfold
~/transcpp-master/scramble
```

