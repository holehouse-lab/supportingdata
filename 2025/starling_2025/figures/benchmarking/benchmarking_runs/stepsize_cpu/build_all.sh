#!/bin/zsh
for i in 1 2 3 4 5 10 15 20 25 30 40 50
do
	starling-benchmark --single-run 200 --batch-size 200 --steps $i --device cpu
	sleep 2
done
	
