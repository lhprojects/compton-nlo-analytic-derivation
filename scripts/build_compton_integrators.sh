#!/usr/bin/env bash
# Build numerical integration engines in a user cache, never in chapter sources.
set -euo pipefail
source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
task_cache="${XDG_CACHE_HOME:-$HOME/.cache}/qftbook/compton-integration"
mkdir -p "$task_cache"
if [[ ! -f "$task_cache/Cuba-4.2.2/libcuba.a" ]]; then
  curl --fail --location https://feynarts.de/cuba/Cuba-4.2.2.tar.gz -o "$task_cache/cuba.tar.gz"
  (cd "$task_cache" && echo '8d9f532fd2b9561da2272c156ef7be5f3960953e4519c638759f1b52fe03ed52  cuba.tar.gz' | sha256sum --check)
  tar -xzf "$task_cache/cuba.tar.gz" -C "$task_cache"
  (cd "$task_cache/Cuba-4.2.2" && ./configure && make -j 1)
fi
# Cuba's archive recipes must run serially to avoid concurrent ar corruption.
g++ -O3 -march=native -funroll-loops "$source_dir/compton_real_integral.cpp" \
  -I"$task_cache/Cuba-4.2.2" "$task_cache/Cuba-4.2.2/libcuba.a" -lm \
  -o "$task_cache/compton_real_integral"
if [[ "${1:-}" == '--with-qcdloop-check' ]]; then
  if [[ ! -d "$task_cache/qcdloop/.git" ]]; then
    git clone https://github.com/scarrazza/qcdloop.git "$task_cache/qcdloop"
  fi
  git -C "$task_cache/qcdloop" checkout --detach 52d94ce9e628c1b1f9784f66a6ef9de285359ae3
  cmake -S "$task_cache/qcdloop" -B "$task_cache/qcdloop/build" -DCMAKE_BUILD_TYPE=Release -DENABLE_FORTRAN_WRAPPER=OFF
  cmake --build "$task_cache/qcdloop/build" -j 4
  g++ -O3 -fPIC -shared "$source_dir/compton_scalar_bridge.cpp" \
    -I"$task_cache/qcdloop/src" -L"$task_cache/qcdloop/build" \
    -Wl,-rpath,"$task_cache/qcdloop/build" -lqcdloop -lquadmath \
    -o "$task_cache/compton_scalar_check.so"
fi
printf 'Built: %s/compton_real_integral\n' "$task_cache"
