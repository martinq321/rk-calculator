#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v em++ >/dev/null || { echo 'Activate the Emscripten SDK first.' >&2; exit 1; }
em++ -O3 -std=c++17 rk_native.cpp \
  --preload-file kernels_fine.bin --preload-file kernels_mesh.bin \
  -sALLOW_MEMORY_GROWTH=1 -sINVOKE_RUN=0 \
  -sEXPORTED_FUNCTIONS='["_main"]' \
  -sEXPORTED_RUNTIME_METHODS='["callMain"]' \
  -o rk_wasm.js
echo 'Serve this directory over HTTP, then open wasm_demo.html.'
