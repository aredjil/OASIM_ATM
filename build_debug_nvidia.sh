#!/usr/bin/env bash

PROFILING=${1:-OFF}

DEST_DIR=builds/debug_nvidia

CURR_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

rm -rf ${CURR_DIR}/${DEST_DIR}
mkdir -p ${CURR_DIR}/${DEST_DIR}

cd ${CURR_DIR}/${DEST_DIR} && cmake ../.. \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DENABLE_PROFILING=${PROFILING}
