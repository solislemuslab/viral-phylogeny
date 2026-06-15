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
```bash

conda activate phylo
cd ~/DeepLearningClaudia/Phyloformer
```




