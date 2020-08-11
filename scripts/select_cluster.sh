#!/bin/bash
# $1: input file
# $2: cluster id to select
# $3: output file

grep "\^$2" $1 > $3
