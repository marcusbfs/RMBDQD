#!/bin/bash

# Possui CMake como depedência

# Nome do programa
prog="RMBDQD"
build="build"

if [ ! -d "$build" ]; then
    echo "Creating build"
    mkdir "$build"
fi

cd "$build"

cmake .. -DCMAKE_INSTALL_PREFIX=..  && make install #VERBOSE=1

cd ../bin
mv program $prog
