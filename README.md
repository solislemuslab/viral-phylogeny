# viral-phylogeny

## Phyloformer
### Install mamba if you want to use it instead of conda
conda install -n base -c conda-forge mamba

### Clone the phyloformer repo
git clone https://github.com/lucanest/Phyloformer.git && cd Phyloformer

### Create the virtual env and install the phyloformer package inside
conda create -n phylo python=3.9 -c defaults && conda activate phylo
pip install -r requirements.txt

### Generate the distance matrix (Change the code)
python infer_alns.py -o data/testdata/pf_matrices models/pf.ckpt data/testdata/msas

### Generate tree from generated matrix
./bin/bin_macos/fastme -i "data/pf66_matrices/PF00066.phy" -o "data/pf66_matrices/PF00066.nwk" --nni --spr


## MSA Transformer
### Install MsaPhylo
git clone https://github.com/Cassie818/MsaPhylo.git
cd MsaPhylo

### Install prerequisite packages
pip install fair-esm --quiet
pip install transformers --quiet
pip install pysam --quiet
pip install Bio
pip install ete3

### Run the model
python MsaPhylo.py
        --i "./data/Pfam/PF00066.fasta" \
        --name 'PF00066' \
        --o "./PF00066_output_tree" \
        --l 2


## CMAPLE


## PhyloGFN
