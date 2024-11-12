#!/usr/bin/env bash

read -r -d '' usage <<-EOF
usage: $(basename "$0") build_dir install_dir

Download, build and install ADS dependencies (fmt, Galois, Catch2)

Positional arguments:
 build_dir       directory where the sources are downloaded and built
 install_dir     installation prefix
EOF

if [[ $# -ne 2 ]]; then
    echo "${usage}"
    exit 1
fi

set -ex

BUILD_DIR=${1}
INSTALL_DIR=${2}

BUILD_TYPE=Release

if [[ -d "${INSTALL_DIR}" ]]; then
  echo "Installation dir ${INSTALL_DIR} exists, skipping installing dependencies"
  exit
fi

mkdir -p "${BUILD_DIR}"
mkdir -p "${INSTALL_DIR}"

CMAKE_BUILD_PARALLEL_LEVEL=$(nproc)
export CMAKE_BUILD_PARALLEL_LEVEL


cd "${BUILD_DIR}"

# Install fmt
FMT_VER=8.0.1
curl -sL https://github.com/fmtlib/fmt/archive/refs/tags/${FMT_VER}.tar.gz -o fmt.tar.gz
tar xzf fmt.tar.gz
rm fmt.tar.gz
mv fmt-${FMT_VER} fmt

mkdir -p fmt/build

cmake \
  -S fmt \
  -B fmt/build \
  -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -D CMAKE_PREFIX_PATH="${INSTALL_DIR}" \
  -D FMT_TEST=OFF \
  -D FMT_DOC=OFF

cmake --build fmt/build
cmake --install fmt/build

# Install Lyra
LYRA_VER=1.5.1
curl -sL https://github.com/bfgroup/Lyra/archive/refs/tags/${LYRA_VER}.tar.gz -o lyra.tar.gz
tar xzf lyra.tar.gz
rm lyra.tar.gz
mv Lyra-${LYRA_VER} Lyra

mkdir -p Lyra/build

cmake \
  -S Lyra \
  -B Lyra/build \
  -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -D CMAKE_PREFIX_PATH="${INSTALL_DIR}"

cmake --build Lyra/build
cmake --install Lyra/build

# Install Galois
GALOIS_VER=6.0
git clone --branch release-${GALOIS_VER} --depth=1 --quiet https://github.com/IntelligentSoftwareSystems/Galois
sed -i '36i #include <optional>' Galois/tools/graph-convert/graph-convert.cpp

mkdir -p Galois/build

cmake \
  -S Galois \
  -B Galois/build \
  -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -D CMAKE_PREFIX_PATH="${INSTALL_DIR}" \
  -D GALOIS_ENABLE_DIST=OFF \
  -D BUILD_TESTING=OFF

cmake --build Galois/build --target galois_shmem
cmake --build Galois/build --target graph-convert
cmake --build Galois/build --target graph-convert-huge
cmake --install Galois/build --component dev
cmake --install Galois/build --component lib
cmake --install Galois/build --component tools

# Install Catch2
CATCH2_VER=3.3.2
curl -sL https://github.com/catchorg/Catch2/archive/refs/tags/v${CATCH2_VER}.tar.gz -o catch2.tar.gz
tar xzf catch2.tar.gz
mv Catch2-${CATCH2_VER} Catch2

mkdir -p Catch2/build

cmake \
  -S Catch2 \
  -B Catch2/build \
  -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -D CMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -D CMAKE_PREFIX_PATH="${INSTALL_DIR}" \
  -D BUILD_TESTING=OFF

cmake --build Catch2/build
cmake --install Catch2/build
