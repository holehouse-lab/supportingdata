# Memory profiling
To profile GPU memory usage, we used `nvidia-smi` to query memory being used. This directory has two files:


* `memory_profile.txt` - this is our nodes captured from running `nvidia-smi` as a function of batch size. We profiled the UNET version of STARLING (as provided in the preprint) and the ViT version of STARLING (as provided in V2 and on GitHub, on both an A4000 and a GTX-1660 (only for the UNET version).

* `memory_prof.csv` - numerical data ONLY for the ViT STARLING implementation which is then plotted in the paper.