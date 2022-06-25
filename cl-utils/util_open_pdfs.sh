#!/bin/bash

cd $1

for f in *.pdf
do
    evince $f
done
