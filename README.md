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

# Run on Linux virtual machine to reduce triangle inequality violation (Even worse: Explained variance is lower)
./bin/bin_linux/fastme \
  -q \
  -i data/testdata/pf_matrices/PF00066.phy \
  -o data/testdata/pf_matrices/PF00066_q.nwk \
  --nni --spr

# Run on MacOS
./bin/bin_macos/fastme \
  -i data/testdata/pf_matrices/PF00066.phy \
  -o data/testdata/pf_matrices/PF00066.nwk \
  --nni --spr
```

Phyloformer is completed!


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
Setup
```bash
# Clone the CMAPLE
git clone https://github.com/iqtree/cmaple.git

# Copy PF00066.fasta to the folder "CMAPLE"
cp MsaPhylo/data/Pfam/PF00066.fasta CMAPLE/
```

Download cmaple-1.1.0-Linux-intel.tar.gz from the url https://github.com/iqtree/cmaple/releases
```bash
# Move the folder to the same level as PF00066.fasta
mv ~/cmaple-1.1.0-Linux-intel.tar.gz ~/DeepLearningClaudia/CMAPLE/
cd DeepLearningClaudia/CMAPLE/

# Extract the folder
tar -xvzf cmaple-1.1.0-Linux-intel.tar.gz
```

```bash
# Run the PF00066.fasta
cd CMAPLE/
./cmaple-1.1.0-Linux-intel/bin/cmaple -aln PF00066.fasta
```
The results of running the commands above:
![EC7892E4-3EC1-43C7-BBBA-BF69120E7601](https://github.com/user-attachments/assets/8f1d776e-3907-46ee-ac1b-869b8031b5c2)

CMAPLE is completed!


## NeuralNJ
Setup
```bash
# Clone the NeuralNJ
git clone https://github.com/DingShizhe/NeuralNJ.git

# Copy PF00066.fasta to the folder "NeuralNJ"
cp MsaPhylo/data/Pfam/PF00066.fasta NeuralNJ/
```

Comment out - phyloformer==0.0.1a4 in environment.yaml
```bash
# Setup the environment
conda env create -f environment.yaml
conda activate NeuralNJ
pip install biopython

# Install raxmlpy separately
pip install pybind11
cd RAxMLpy
python setup.py build
python setup.py install --user

export LD_LIBRARY_PATH="$(pwd)/build_plllib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$(pwd)/build_raxmllib:$LD_LIBRARY_PATH"
```

```bash
# Create a new python file "convert_to_phylip.py" to convert fasta file to phy file
python convert_to_phylip.py
mv PF00066.phy examples/
```

Modify config file: Replace the previous instance_path with 
```bash
instance_path: "./examples/PF00066/"
```

```bash
# Run the PF00066.fasta
python finetune_rl_search.py \
  --config ./config/finetune_reinforce_search_example.yaml \
  --infer_opt Argmax
```


## PhyloGFN
Setup
```bash
# Clone the PhyloGFN
git clone https://github.com/zmy1116/phylogfn.git

# Copy PF00066.fasta to the folder "PhyloGFN"
cp MsaPhylo/data/Pfam/PF00066.fasta PhyloGFN/

# Setup the environment
conda create -n phylogfn python=3.10
conda activate phylogfn
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
conda install anaconda::docopt
conda install etetoolkit::ete3
conda install matplotlib tqdm dill fvcore iopath docopt
```


