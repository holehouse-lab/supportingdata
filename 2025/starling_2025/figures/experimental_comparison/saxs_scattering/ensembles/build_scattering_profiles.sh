#!/bin/zsh
source ~/.zshrc

for file in *.starling; do
    name="${file%.starling}"
    echo "$name"
    if [ ! -d $name ]
    then
	mkdir $name
    
	cd $name
	mv ../${name}.starling .
	mv ../${name}_STARLING.xtc .
	mv ../${name}_STARLING.pdb .
	
	# NB - autoSCTR is a holehouselab tool that wraps around FOXS for calculating
	# average scattering curves.
	autoSCTR --pdb ${name}_STARLING.pdb --xtc ${name}_STARLING.xtc -n 400
	cd ..
    else
	echo "Already found ensemble for $name"
    fi
    
done
