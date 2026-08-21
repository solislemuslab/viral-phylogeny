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

## IQ-Tree
```bash
# Setup the environment
cd IQTREE_results
module load Miniforge3-26.1.1-3 # Load for conda
conda create -n iq_tree python=3.9 -c defaults
conda activate iq_tree
conda install -c bioconda iqtree

cp ../OB45_concat/OB45_supermatrix.fasta .
iqtree -s OB45_supermatrix.fasta -m GTR+G -nt AUTO
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

### Change from CS department virtual machine to discovery building virtual machine
Setup the environment
```bash
conda activate phylo
cd /mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/Phyloformer

mkdir data/testdata/ob45_msas/
cp ../OB45_concat/OB45_supermatrix.fasta data/testdata/ob45_msas/
```

Generate the distance matrix
```bash
python infer_alns.py \
  -o data/testdata/ob45_concat_matrices \
  models/pf.ckpt \
  data/testdata/ob45_msas
```

Generate tree using the distance matrix above
```bash
# Run on Linux virtual machine
./bin/bin_linux/fastme \
  -i data/testdata/ob45_concat_matrices/OB45_supermatrix.phy \
  -o data/testdata/ob45_concat_matrices/OB45_supermatrix.nwk \
  --nni --spr
```
Done


## MSA Transformer
### Embedding tree
Setup the environment
```bash
conda activate msaphylo
cd ~/DeepLearningClaudia/MsaPhylo

mkdir data/OB45_msas
mkdir OB45_embedding_trees

cp ../OB45_concat/OB45_supermatrix.fasta data/OB45_msas/
```

Generate the trees
```bash
python MsaPhylo.py \
  --i "./data/OB45_msas/OB45_supermatrix.fasta" \
  --name OB45 \
  --o "./OB45_output_tree" \
  --l 2
```

Failure:
```bash
Downloading: "https://dl.fbaipublicfiles.com/fair-esm/models/esm_msa1b_t12_100M_UR50S.pt" to /mnt/ws/home/jgao/.cache/torch/hub/checkpoints/esm_msa1b_t12_100M_UR50S.pt
Downloading: "https://dl.fbaipublicfiles.com/fair-esm/regression/esm_msa1b_t12_100M_UR50S-contact-regression.pt" to /mnt/ws/home/jgao/.cache/torch/hub/checkpoints/esm_msa1b_t12_100M_UR50S-contact-regression.pt
Traceback (most recent call last):
  File "/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/MsaPhylo/MsaPhylo.py", line 120, in <module>
    main()
  File "/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/MsaPhylo/MsaPhylo.py", line 116, in main
    embtree.build_emb_tree()
  File "/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/MsaPhylo/MsaPhylo.py", line 61, in build_emb_tree
    plm_embedding = self.get_embedding()
  File "/mnt/dv/wid/projects4/SolisLemus-viral-phylogeny/DeepLearningClaudia/MsaPhylo/MsaPhylo.py", line 46, in get_embedding
    raise ValueError("It exceeds the capacity of the MSA transformer!")
ValueError: It exceeds the capacity of the MSA transformer!
```

### Attention tree
Create a new python file "run_ob45_attention_trees.py"
```bash
cd data/Pfam/
rm -rf *.fasta
cd ..
cd ..
cp data/OB45/*.fasta data/Pfam/

python run_ob45_attention_trees.py
```


