#!/bin/bash

for file in `ls examples`
do
    echo " "
    echo $file
    ./solve examples/$file
done