# PhyloBench data

```bash
# Upload the dataset "OB45.tar.gz" to the virtual machine
scp ./DeepLearningClaudia/OB45.tar.gz jiayig@best-linux.cs.wisc.edu:/home/jiayig/DeepLearningClaudia/
mkdir OB45
tar -xvzf OB45.tar.gz -C OB45

# Setup the environment
conda create -n phylobench_ob45 python=3.10 -y
conda activate phylobench_ob45
pip install biopython pandas openpyxl tqdm
```

Create a python file called preprocess_ob45.py
```bash
# Run the python file
python preprocess_ob45.py
```


## Phyloformer
Setup the environment
```bash
conda activate phylo
cd ~/DeepLearningClaudia/Phyloformer

mkdir data/testdata/ob45_msas/
cp ~/DeepLearningClaudia/OB45_concat/OB45_supermatrix.fasta data/testdata/ob45_msas/
```

Generate the distance matrix (Error: memory allocation problem)
```bash
python infer_alns.py \
  -o data/testdata/ob45_concat_matrices \
  models/pf.ckpt \
  data/testdata/ob45_msas
```

First try only one fasta file (This works successfully)
```bash
mkdir data/testdata/ob45_one_msa
mkdir data/testdata/ob45_one_matrix
mkdir data/testdata/ob45_one_tree
cp /home/jiayig/DeepLearningClaudia/OB45_cleaned/OB_PF01479_1.fasta data/testdata/ob45_one_msa/

# Generate the distance matrix 
python infer_alns.py \
  -o data/testdata/ob45_one_matrix \
  models/pf.ckpt \
  data/testdata/ob45_one_msa

# Generate tree using the distance matrix above
./bin/bin_linux/fastme \
  -i data/testdata/ob45_one_matrix/OB_PF01479_1.phy \
  -o data/testdata/ob45_one_tree/OB_PF01479_1.nwk \
  --nni --spr
```

SO change to use 603 seperate files rather than the concatenated OB45 dataset
```bash
conda activate phylo
cd ~/DeepLearningClaudia/Phyloformer

rm -rf data/testdata/ob45_msas
rm -rf data/testdata/ob45_matrices
rm -rf data/testdata/ob45_trees

mkdir data/testdata/ob45_msas
mkdir data/testdata/ob45_matrices
mkdir data/testdata/ob45_trees

cp ~/DeepLearningClaudia/OB45_cleaned/*.fasta data/testdata/ob45_msas/
```

Generate the distance matrix
```bash
python infer_alns.py \
  -o data/testdata/ob45_matrices \
  models/pf.ckpt \
  data/testdata/ob45_msas
```

Generate tree using the distance matrix above
```bash
for phy in data/testdata/ob45_matrices/*.phy; do
    base=$(basename "$phy" .phy)
    echo "Running FastME for $base"

    ./bin/bin_linux/fastme \
      -i "$phy" \
      -o "data/testdata/ob45_trees/${base}.nwk" \
      --nni --spr
done
```


## MSA Transformer
### Embedding tree
Setup the environment
```bash
conda activate msaphylo
cd ~/DeepLearningClaudia/MsaPhylo
```



