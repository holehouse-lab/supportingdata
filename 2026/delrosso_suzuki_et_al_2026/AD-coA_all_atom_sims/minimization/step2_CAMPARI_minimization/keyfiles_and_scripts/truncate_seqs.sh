#!/usr/bin/zsh

### ARGUMENTS ###
nreps=3
input_pdb_name=gromacs_minim.pdb
output_pdb_name=truncated.pdb
campari_minimize_keyfile=minimize.key

#################

currdir=$(pwd)

rm run_min_campari_parallel.sh
touch run_min_campari_parallel.sh

while read -r line; do
	echo $line
	name=$(echo $line | awk '{print $1}')
	motif_seq=$(echo $line | awk '{print $2}')
	motif_start=$(echo $line | awk '{print $3}')
	motif_end=$(echo $line | awk '{print $4}')
	idr_chain=$(echo $line | awk '{print $5}')
	fd_seq=$(echo $line | awk '{print $6}')

	# Create seq.in file
	echo "${name} ${motif_seq}" > seq.fasta

	if [ "$idr_chain" = "0" ]; then
		python make_seq_in_file.py $fd_seq seq.fasta seq.in --idr-caps --idr-first
	else
		python make_seq_in_file.py $fd_seq seq.fasta seq.in --idr-caps
	fi

	for i in {1..$nreps}; do
		cd $i
		cp ../seq.in .
		rm campari_minimized* ENERGY.dat POLYMER.dat minimize.key
		echo "cd ${currdir}/${i}; cp ../${campari_minimize_keyfile} .; campari3 -k ${campari_minimize_keyfile}" >> ../run_min_campari_parallel.sh

		if [ "$idr_chain" = "0" ]; then
			python modify_complex_pdb.py --pdb $input_pdb_name --out $output_pdb_name --clip-idx $motif_start $motif_end --idr-first
		else
			python modify_complex_pdb.py --pdb $input_pdb_name --out $output_pdb_name --clip-idx $motif_start $motif_end
		fi

		cd ..
	done

done < sequence_data.txt
