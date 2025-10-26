#!/usr/bin/env bash


DEST_DIR=builds/debug_gcc

CURR_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

rm -rf ${CURR_DIR}/${DEST_DIR}
mkdir -p ${CURR_DIR}/${DEST_DIR}

cd ${CURR_DIR}/${DEST_DIR} && cmake -D CMAKE_BUILD_TYPE=Debug -DENABLE_PROFILING=ON -DCMAKE_PREFIX_PATH=/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/netcdf-fortran-4.6.1-nvymynen3aoyohew4d7t646sjig4ufl6 ../..

