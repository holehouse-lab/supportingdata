#!/bin/zsh

##
## ZENODO assembly script
## J. Alston & ASH
##
## Last updated 2022-02-10
##
## Changelog
## V1 - initila version
##
##

## ABOUT
# This script constructs the Zenono-shared .zip file with our trajectories for the manuscript
#
# The disordered N-terminal tail of SARS CoV-2 Nucleocapsid protein forms a dynamic complex with RNA
# Jasmine Cubuk, Jhullian J. Alston, J. Jeremías Incicco, Alex S. Holehouse, Kathleen B Hall, Melissa Brereton, Andrea Soranno
#
# For questions pertaining to simulations please contact alex.holehouse@wustl.edu or Jhulian Alston
#
# What is this script doing?
#
# This script was used to assemble the final shared set of trajectory data we provide accompanying the manuscript above. Note we have 2 more
# independent replicas (80,000 frames each) if you want more data, but the sampling and volume of data here reveals essentially identical results
# to the average over the three datasets so, in the interest of keeping things tractable we are only sharing 1/3 of the data. If you want the
# remaining 2/3 please let us know...
#
# We provide this script to show EXACTLY how the data were assembled from within the Holehouse lab's filesystem, making it easy to understand how
# and where simulations were run.


## A request
#
# We are still working on a follow-up manuscript using these data, so we kind ask IF you are interested in re-analyzing this and publishing your
# own paper then (1) That's GREAT but (2) Please consider shooting Alex or Jhullian an email just so we can coordinate and make sure we're not
# doubling up on the same question(s). This would be greatly appreciated.
#
#


##
## FIRST: NTD, RBD, and NTD-RBD sims
##

for i in NTD RBD NTDRBD
do
    for u in 0  10  12  15  17  180  20  25  30  35  40
    do

	# make the output directory (-p flag means full path is constructed if any of
	# the subdirectories don't YET exist)
	mkdir -p ${i}/${u}/1/ 

	# copy topology file 
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/equilibrated.pdb ${i}/${u}/1/

	# copy trajectory file (equilibrated just means for 0.2% has been discarded as
	# equilibration as per methods
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/equilibrated.xtc ${i}/${u}/1/


	## NOTE the three operations below do error for RNA len = 0 but that's OK. we have
	# not set 'set -e' so shouldn't be an issue
	
	# copy the center-of-mass distances. Note in our notebooks we actually calculate
	# these directly from the sims, but given the info needed for computing apparent
	# K_A is that COM distance we provided this info directly
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/com_distance.csv ${i}/${u}/1/

	# copy the bound trajectories (note the '5' here reflects the fact we required 5 consecutive
	# frames to count something as bound. Truthfully, the binding affinities are not particularly
	# sensitive to this number but you can play around for yourself with the code if that's your jam
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/bound_5.xtc ${i}/${u}/1/bound.xtc

	# copy the unbound trajectories, same logic as operation above
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/unbound_5.xtc ${i}/${u}/1/unbound.xtc


    done
done

## Next do omicron variants with rU25

for i in Omicron P13L D3133
do
    for u in 25
    do

	# make the output directory (-p flag means full path is constructed if any of
	# the subdirectories don't YET exist)
	mkdir -p ${i}/${u}/1/ 

	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/equilibrated.pdb ${i}/${u}/1/
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/equilibrated.xtc ${i}/${u}/1/
	cp /work/alstonj/2022/Binding/FINAL/${i}/${u}/1/com_distance.csv ${i}/${u}/1/

    done
done

echo "Directories complete"
	
