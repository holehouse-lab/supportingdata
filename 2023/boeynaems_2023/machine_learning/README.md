# To replicate the analysis please run the following commands:
conda create -n venom python=3.11 -y
conda activate venom

pip install jax-unirep seaborn matplotlib protfasta

python train_model.py
