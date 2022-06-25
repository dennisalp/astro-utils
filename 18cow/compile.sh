#!/bin/bash

if [ $1 == "debug" ]
then
    CC_="mpicc"
    FLAGS="-Werror -Wall -Wextra -g3 -O0 -DDEBUG"
    SANITIZER="-fno-omit-frame-pointer -fsanitize=address,undefined,float-divide-by-zero,float-cast-overflow,alignment,nonnull-attribute,null,signed-integer-overflow"
    mkdir -p bin out
elif [ $1 == "release" ]
then
    CC_="mpicc"
    FLAGS="-Ofast -march=native"
    mkdir -p bin out
elif [ $1 == "tegner" ]
then
    CC_="mpiicc"
    FLAGS="-Ofast -march=native"
    mkdir -p bin out
elif [ $1 == "clean" ]
then
    rm -rf bin/
    exit 0
else
    echo "Invalid argument"; exit 1
fi

$CC_ $FLAGS $SANITIZER src/cow.c -o bin/cow

#CFLAGS=-I/Users/silver/mpi/include
#EXTRA_FLAGS=-Werror -Wall -Wextra -g3 -O0 -DDEBUG
#SANITIZER_FLAGS=-fno-omit-frame-pointer -fsanitize=address,undefined,float-divide-by-zero,float-cast-overflow,alignment,nonnull-attribute,null,signed-integer-overflow
#INCLUDES=-Isrc/include
#DEPENDENCY_FLAGS=-MMD -MP
#
#CFLAGS=-I/Users/silver/mpi/include
#EXTRA_FLAGS=-Ofast -march=native
#SANITIZER_FLAGS=
#INCLUDES=-Isrc/include
#DEPENDENCY_FLAGS=-MMD -MP
