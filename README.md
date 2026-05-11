# Viral-phylogeny

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
The file `PF00066.fasta` is copied into the `data/testdata/msas`, and then generate the distance matrix
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
## this generated a warning about violation to triangle inequality

# Run on Linux virtual machine to reduce triangle inequality violation (Even worse: Explained variance is lower)
./bin/bin_linux/fastme \
  -q \
  -i data/testdata/pf_matrices/PF00066.phy \
  -o data/testdata/pf_matrices/PF00066_q.nwk \
  --nni --spr
## same warning!

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
Create a new python file "run_attention_pf00066.py" for attention tree. Name of output file is inside this python script.
This script is in the `scripts` folder.

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
Note this requires compilation. Note that there is a compilation step for DNA and another one for AA.

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
```

Add the code below to include_dirs=[...] in setup.py
```bash
cd RAxMLpy
"raxml-ng/libs/coraxlib/src",
```

```bash
# Install raxmlpy separately
pip install pybind11
python setup.py build
python setup.py install

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

## Phyloformer2
Setup
```bash
# Clone the Phyloformer2
git clone --recursive git@gitlab.in2p3.fr:deelogeny/wp1/phyloformer-2.git

# Copy PF00066.fasta to the folder "Phyloformer2/"
cp MsaPhylo/data/Pfam/PF00066.fasta Phyloformer2/

# Setup the environment with conda
conda create -n phylo2 python=3.10
conda activate phylo2
cd /home/jiayig/DeepLearningClaudia/Phyloformer2
conda install --yes --file requirements.txt
python path/to/pf2/script.py

# Setup the environment with uv 
cd path/to/pf2/repo
uv init
uv add -r requirements.txt
uv run path/to/pf2/script.py
```



# COVID data

In GISAID website, you can upload the `msaa314_Supplementary_Data.xls` file that we get from the paper [Supp Mat](https://academic.oup.com/mbe/article/38/5/1777/6030946#304310196) and you download two fasta files (because you cannot download more than 10,000 sequences at once). However, we did not follow this approach.

Jiayang followed the information in the paper on how data was collected:
```
We downloaded the raw data from gisaid.org on May 5, 2020. It contained 16,453 full-genome (⁠>29,000bp) raw sequences with high coverage. High-coverage sequences are defined by GISAID as sequences containing less than 1% Ns (undetermined characters), less than 0.05% unique amino acid mutations, and no insertions/deletions unless these have been verified by the submitter.
```

There are two fasta files:
- `gisaid_hcov-19_after_2020_03_30.fasta`
- `gisaid_hcov-19_to_2020_03_30.fasta`

Four outgroups are included after (see folder `outgroup` in the google drive).

**Note:** Joyce extracted 300 sequences for testing, but we want to expand to all. And Jiayang thinks we should download the data better.

```bash
# Setup the environment
conda create -n covid_data python=3.10 -y
conda activate covid_data
pip install biopython pandas openpyxl tqdm

# Preprocess the COVID dataset
cd COVID_data_fasta/
cat gisaid_hcov-19_to_2020_03_30.fasta \
    gisaid_hcov-19_after_2020_03_30.fasta \
    > covid_all.fasta
```

Create remove_duplicates.py to remove duplicated in the new-created fasta file
```bash
python remove_duplicates.py
```

Create select_covid_300_stratified.py to select 300 sequences to run for different models
```bash
python select_covid_300_stratified.py
```
Below is the distribution of 300 sequences

Selected counts by continent:

continent

- Europe           146
- North America     93
- Oceania           30
- Asia              27
- Africa             2
- South America      2

Add outgroups to the covid_300_stratified.fasta
```bash
cat covid_300_stratified.fasta outgroup/*.fasta > covid_300_stratified_with_outgroup.fasta
```

Aligned the 300 sequences
```bash
conda install -c bioconda mafft
mafft --auto covid_300_stratified_with_outgroup.fasta > covid_300_aligned.fasta
```

### IQ-Tree
```bash
# Run the IQ-Tree
cd COVID_data_fasta
conda activate covid_data
conda install -c bioconda iqtree
iqtree -s covid_300_aligned.fasta -m GTR+G -nt AUTO

# Move the results of IQ-Tree to a separate folder
mkdir IQTREE_results
mv COVID_data_fasta/covid_300_aligned.fasta.* IQTREE_results/ # covid_300_aligned.fasta.treefile is the final tree
```

### CMAPLE
```bash
# Run the CMAPLE
conda deactivate
cd CMAPLE/
./cmaple-1.1.0-Linux-intel/bin/cmaple -aln ~/DeepLearningClaudia/COVID_data_fasta/covid_300_aligned.fasta

# Move the results of CMAPLE to a separate folder
mkdir CMAPLE_results
cd ..
mv COVID_data_fasta/covid_300_aligned.fasta.* CMAPLE/CMAPLE_results/ # covid_300_aligned.fasta.treefile is the final tree
```

### NeuralNJ
```bash
# Run the NeuralNJ
conda deactivate
conda activate NeuralNJ
cd NeuralNJ/
```

Change the name of file in convert_to_phylip.py to convert fasta file to phy file
```bash
cp ~/DeepLearningClaudia/COVID_data_fasta/covid_300_aligned.fasta .
python convert_to_phylip.py
mv covid_300_aligned.phy examples/
```

Modify config file: Replace the previous instance_path with
```bash
instance_path: "./examples/covid_300_aligned/"
```

```bash
# Run the PF00066.fasta
python finetune_rl_search.py \
  --config ./config/finetune_reinforce_search_example.yaml \
  --infer_opt Argmax
```

```bash
# Move the results of NeuralNJ to a separate folder
mkdir NeuralNJ_results
cd ..
mv COVID_data_fasta/covid_300_aligned.fasta.* NeuralNJ/NeuralNJ_results/ # covid_300_aligned.fasta.treefile is the final tree
```

### Phyloformer
```bash
# Run the phyloformer
conda deactivate
conda activate phylo
cd Phyloformer/
mkdir data/testdata/covid_msas
cp ~/DeepLearningClaudia/COVID_data_fasta/covid_300_aligned.fasta data/testdata/covid_msas/
python infer_alns.py \
  -o data/testdata/covid_matrices \
  models/pf.ckpt \
  data/testdata/covid_msas

# Move the results of phyloformer to a separate folder
mkdir Phyloformer_results
cd ..
mv COVID_data_fasta/covid_300_aligned.fasta.* Phyloformer/Phyloformer_results/ # covid_300_aligned.fasta.treefile is the final tree
```

## Visualization
Setup the environment
```bash
conda create -n visualization python=3.9 -c defaults
conda install -c conda-forge ete3 biopython matplotlib
```

Change the path in python file compare.py to compare IQ-Tree and CMAPLE
```bash
python compare.py
```

Below are the results:
RF distance: 436
Maximum RF distance: 602
Normalized RF distance: 0.7242524916943521
Number of common leaves: 304
Common taxa: 304





