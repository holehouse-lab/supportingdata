# Instructions for running comparison with IDPForge
##### Last updated 2025-06-03

The documentation below outlines the steps required to get the IDPForge Docker container up and running on a Linux workstation with NVIDIA GPU drivers installed and Docker configured to utilize those drivers.

### Testing GPUs are accessible

As of June 2025, the following command works for testing if GPUs are accessible:

	docker run --rm --gpus all nvidia/cuda:12.2.0-base-ubuntu22.04 nvidia-smi
	
This should return the `nvidia-smi` output if all is well. This is the first important thing to check.

### Getting the IDPForge container built
As described on the IDPForge documentation, you can build the container by cloning the idpforge repo and running:

	docker build -t idpforge .

NOTE: This takes 3-4 hours to build (most of which is step 8/15).

### Fixing the IDPForge container
There are a number of small issues that prevent the IDPForge container from running out of the box. Below, we document the changes we made to get it running.

#### 1. Add explicit definitions for ESMFold StructureModule initialization

When trying to run in the container, the first error to pop up is:

	File "/opt/idpforge/sample_idp.py", line 114, in <module>
		main(args.ckpt_path, args.output_dir, args.sample_cfg,
	File "/opt/idpforge/sample_idp.py", line 24, in main
		model = IDPForge(settings["diffuse"]["n_tsteps"],
	File "/opt/idpforge/idpforge/model.py", line 69, in __init__
		self.trunk = FoldingTrunk(**cfg.trunk)
	File "/opt/idpforge/esm/esmfold/trunk.py", line 150, in __init__
		self.structure_module = StructureModule(**self.cfg.structure_module)
	TypeError: __init__() missing 3 required positional arguments: 'trans_scale_factor', 'epsilon', and 'inf'
	
This is because the StructureModule object – part of the ESMFold architecture – initializes based on the configuration data passed through the `sample.yaml` file, and these three fields are missing from that file. To fix this, we need to edit the sample file to include those fields. 

Specifically, we must change the sample file being used (e.g., `configs/sample.yml`) from:


	 trunk:
	    num_blocks: 2
	    sequence_state_dim: 128
	    pairwise_state_dim: 64
	    sequence_head_width: 32
	    pairwise_head_width: 32
	    max_recycles: 3
	    recycle_min_bin: 3.375
	    recycle_max_bin: 39.375
	    structure_module:
	      c_s: 256
	      c_z: 64
	      c_ipa: 16
	      c_resnet: 128
	      no_heads_ipa: 8
	      no_qk_points: 4
	      no_v_points: 8
	      dropout_rate: 0.1
	      no_blocks: 4
	      no_transition_layers: 1
	      no_resnet_blocks: 2
	      no_angles: 7

To (see added bottom three lines)

	 trunk:
	    num_blocks: 2
	    sequence_state_dim: 128
	    pairwise_state_dim: 64
	    sequence_head_width: 32
	    pairwise_head_width: 32
	    max_recycles: 3
	    recycle_min_bin: 3.375
	    recycle_max_bin: 39.375
	    structure_module:
	      c_s: 256
	      c_z: 64
	      c_ipa: 16
	      c_resnet: 128
	      no_heads_ipa: 8
	      no_qk_points: 4
	      no_v_points: 8
	      dropout_rate: 0.1
	      no_blocks: 4
	      no_transition_layers: 1
	      no_resnet_blocks: 2
	      no_angles: 7
	      trans_scale_factor: 10
	      epsilon: 0.00000001
	      inf: 100000000.0


#### 2. Fix the model keys 
In `sample_idp.py`, we need to update the model reading code because the model keys in the provided weights file don't match the expected model keys. Specifically, we read the model info by adding the following code:


	fixed_sd = {}
	for k, v in sd["ema"]["params"].items():
	    if "linear_q_points.weight" in k:
	        k = k.replace("linear_q_points.weight", "linear_q_points.linear.weight")
	    elif "linear_q_points.bias" in k:
	        k = k.replace("linear_q_points.bias", "linear_q_points.linear.bias")
	    elif "linear_kv_points.weight" in k:
	        k = k.replace("linear_kv_points.weight", "linear_kv_points.linear.weight")
	    elif "linear_kv_points.bias" in k:
	        k = k.replace("linear_kv_points.bias", "linear_kv_points.linear.bias")
	    fixed_sd[k] = v
	
	model.load_state_dict(fixed_sd)

#### 3. Replace `example_data.pkl` with a real secondary structure file
The `sample.yml` file, which defines the settings used for prediction, references a file called `data/example_data.pkl`. This file is unfortunately not in the `/data` directory, but can be found in the [figshare repo](https://doi.org/10.6084/m9.figshare.28414937). However, this is a pickle file, and the code reading here expects a text file; specifically, reading this file raises:

	Traceback (most recent call last):
	  File "/opt/idpforge/sample_idp.py", line 112, in <module>
	    main(args.ckpt_path, args.output_dir, args.sample_cfg,
	  File "/opt/idpforge/sample_idp.py", line 80, in main
	    ss = f.read().split("\n")
	  File "/opt/conda/lib/python3.9/codecs.py", line 322, in decode
	    (result, consumed) = self._buffer_decode(data, self.errors, final)
	UnicodeDecodeError: 'utf-8' codec can't decode byte 0x80 in position 0: invalid start byte

This is indicative of the fact that the code is reading, expecting a text stream, but gets a binary (pickle) file.

We therefore have to update the `sample.yml` file to read instead a text secondary structure file, which itself is generated using `prep_sec.py` file. 

The secondary structure file can be generated by moving to the `idpforge/utils` directory and then running the command:

	python prep_sec.py <SEQUENCE> data/example_data.pkl <NUMBER_OF_CONFORMERS> <PATH_OUTPUT_FILE>
	
> NB this is not the command provided in the docs, but the argparser does not define the flagged arguments. 

For example, to generate a secondary structure file for 100 conformers A1_LCD_10R sequence:

	python prep_sec.py GSMASASSSQGGSSGSGNFGGGGGGGFGGNDNFGGGGNFSGSGGFGGSGGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGSSSGPYGGGGQYFAKPGNQGGYGGSSSSSSYGSGGGF ../../data/example_data.pkl 100 ../../data/A1_LCD_10R_sample_100.txt
	
We must generate one of these for each sequence being predicted, and then make sure `sample.yml` is updated to point at this file.

#### 4. Install missing/outdated software
We also need to install some missing dependencies inside the docker container

	# needed for structure relaxing
	conda install -c conda-forge openmm
	
	# needed because the original pdbfixer is too old (requires openmm < 7.6 
	conda install -c conda-forge pdbfixer

#### 5. Fix the broken PDBWriter in OpenFold
OpenFold's pdbwriter appears to be broken; specifically, in the `openfold.np.protein.py` file, if chain_index is set to None in the `to_pdb()` function this raises an exception:

	AttributeError: 'NoneType' object has no attribute 'astype'

I honestly don't know how this ever could have worked (note `chain_index` is an optional argument which defaults to `None`), but we can fix it by introducing the following change into the actual raw Python code in the container. i.e., we edit:

	/opt/conda/lib/python3.9/site-packages/openfold-2.2.0-py3.9-linux-x86_64.egg/openfold/np/protein.py
	
Which obviously you should never do in real life, but this is container life, so... 	
Specifically, we add the following code to where chain_index gets initialized in the `to_pdb()` function (around line 342):

	chain_index = prot.chain_index
	
	if chain_index is None:
    	chain_index = np.zeros_like(prot.residue_index, dtype=np.int32)
	else:
    	chain_index = chain_index.astype(np.int32)

#### 6. Get the weights
Finally, we need to get the final weights from the[figshare repo](https://doi.org/10.6084/m9.figshare.28414937), specifically the file called `mdl.ckpt`. I copied these into the container into a directory I created, called `/weights` inside `/opt/idpforge`.					
		
## Test-driving 		
After completing these steps, the container is ready to run. To test this out, I suggest first starting the container in an interactive way that mounts your current directory directly as a mountpoint:

	docker run -it --rm --gpus all -v "$(pwd):/mount" idpforge_ash:latest
	
This will drop you into the `/opt/idpforge` directory, and from there you can check things are working by:

1. Generating the secondary structure file for your sequence of interest.
2. Updated `sample.yaml	` to point to the correct secondary structure (ss) file.
3. Running the actual IDPforge sample command.

Points 1 and 2 are described above, so to run the command, the following should work.

	python sample_idp.py weights/mdl.ckpt test <sample_file_here> --nconf 100 --cuda
	
Doing this, I get 10-20 seconds per conformer for ~130 residue sequence (note relaxation is only part of the prediction time) to generate fully relaxed conformations on an A4500 Nvidia GPU, which seems comparable to the numbers reported by the authors.	

	
## predict_all.sh
The `predict_all.sh` script is then placed inside the idpforge directory (`/opt/idpforge/`), where it can be executed to read the `name_seq.tsv` file, which contains the name and sequence mappings for the 


nohup docker run --rm --gpus all -v "$(pwd):/mount" idpforge_ash:latest bash /opt/idpforge/predict_all.sh &

	# -------------------------------------	
	#!/bin/bash
	
	# set the number of conformations to generate
	N_CONFS=10
	
	# loops over the file at the mountpoint called name_seq.tsv, which expects to be
	# of the format
	#
	# name sequence
	#
	# with a single space separating the name and the sequence
	
	while IFS=$' ' read -r name seq; do
	    echo "On: $name"
	
	    # update the sample file
	    sed -i.bak "s/^sequence: .*/sequence: $seq/" configs/sample2.yml
	
	    # build outdir
	    mkdir -p /mount/predictions/${name}
	
	    # remove any existing pdb files so we generate new predictions
	    rm /mount/predictions/${name}/*_relaxed.pdb
	
	    # get start-time
	    START=$(date +%s)
	
	    # build secondary stucture file
	    python idpforge/utils/prep_sec.py ${seq} data/REAL_example_data.pkl $N_CONFS data/sec_100.txt
	
	    # run prediction
	    python sample_idp.py weights/mdl.ckpt /mount/predictions/${name} configs/sample2.yml --nconf ${N_CONFS} --cuda
	
	    # get end time
	    END=$(date +%s)
	
	    # calculate seconds for end-start
	    DELTA=$((END - START))
	    d=$(date)
	
	    # write performance info
	    echo "$N_CONFS, $DELTA, $d, $name, $seq" >> /mount/predictions/${name}/performance.txt
	
	done < /mount/name_seq.tsv	