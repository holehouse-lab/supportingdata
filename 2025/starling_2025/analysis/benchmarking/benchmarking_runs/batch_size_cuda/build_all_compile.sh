#!/bin/zsh

# script to test with --compile flag set
for i in 10 20 30 40 60 80 100 120 140 160 180 200 220
do
	starling-benchmark --single-run 500 --batch-size $i --compile
done
	
