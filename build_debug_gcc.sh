#!/usr/bin/env bash


DEST_DIR=builds/debug_gcc

CURR_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

rm -rf ${CURR_DIR}/${DEST_DIR}
mkdir -p ${CURR_DIR}/${DEST_DIR}

cd ${CURR_DIR}/${DEST_DIR} && cmake ../.. 
