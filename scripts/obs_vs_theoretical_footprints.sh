#!/bin/bash
# $1: theoretical file
# $2: observed
# $3: output file

awk '{print $4"\ttheoretical"}' $1 > $3 
awk '{print $4"\tobserverd"}' $1 >> $3

