#!/bin/zsh
for i in 50 100 150 200 250 300 350 400 450 500
do
	starling-benchmark --single-run 500 --batch-size $i --device cuda
done
	
