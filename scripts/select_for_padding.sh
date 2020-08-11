#!/bin/bash
awk '{if (NF == pad) {print $0}}' pad=$3 $1 > $2 
