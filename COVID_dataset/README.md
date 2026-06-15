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

### Phyloformer (not run right now because this is for aminoacids)
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
