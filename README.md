# viral-phylogeny

## Phyloformer
Setup
```bash
# Clone the phyloformer repo
git clone https://github.com/lucanest/Phyloformer.git

# Copy PF00066.fasta to the folder "Phyloformer"
cp MsaPhylo/data/Pfam/PF00066.fasta data/testdata/msas/

# Create the virtual env and install the phyloformer package inside
cd Phyloformer
conda create -n phylo python=3.9 -c defaults
conda activate phylo
```

Change the pandas version from 3.5.1 to pandas==2.0.3 in the requirements.txt 
```bash
pip install -r requirements.txt
```

Generate the distance matrix
```bash
python infer_alns.py \
  -o data/testdata/pf_matrices \
  models/pf.ckpt \
  data/testdata/msas
```

Generate tree using the distance matrix above
```bash
# Run on Linux virtual machine
./bin/bin_linux/fastme \
  -i data/testdata/pf_matrices/PF00066.phy \
  -o data/testdata/pf_matrices/PF00066.nwk \
  --nni --spr

# Run on MacOS
./bin/bin_macos/fastme \
  -i data/testdata/pf_matrices/PF00066.phy \
  -o data/testdata/pf_matrices/PF00066.nwk \
  --nni --spr
```


## MSA Transformer
### Embedding tree
Install MsaPhylo
```bash
git clone https://github.com/Cassie818/MsaPhylo.git
cd MsaPhylo
```

Create an environment
```bash
conda create -n msaphylo python=3.10
conda activate msaphylo
```

Install prerequisite packages
```bash
pip install torch torchvision torchaudio
pip install fair-esm --quiet
pip install transformers --quiet
pip install pysam --quiet
pip install Bio
pip install ete3
```

Run the model
```bash
python MsaPhylo.py \
  --i "./data/Pfam/PF00066.fasta" \
  --name PF00066 \
  --o "./PF00066_output_tree" \
  --l 2
```

### Attention tree
Create a new python file "run_attention_pf00066.py" for attention tree
Build the tree
```bash
mkdir -p Attentions/Pfam
python run_attention_pf00066.py
```
The generated trees of using attention tree are stored under the directory "Trees".

MSA Transformer is completed!


## CMAPLE


## PhyloGFN
