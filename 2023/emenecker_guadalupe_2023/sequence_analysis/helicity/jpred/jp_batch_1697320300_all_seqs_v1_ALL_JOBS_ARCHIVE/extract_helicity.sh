#!/bin/zsh

for i in jp_*
do
	cd $i
	cat ${i}.jnet | grep JNETPROPH | awk -F ":" {' print $2 '} > helicity_scores.csv
	cd ..
done