ARStack
========================================

## Environment Setup
The code has been tested successfully using Python 3.8; thus we suggest using this version or a later version of Python. A typical process for installing the package dependencies involves creating a new Python virtual environment.

To install the required packages, run the following:
```console
pip install .
```
The following command downloads the required files, including MSA inputs, Alphafold and Rossetta predicted structures and ground truth structures. 
```console
python download_data.py
```
## running codes
Generate processed MSA input data from the current MSA data in txt file format, which consists of the following steps:
- txt to csv 
```console
python3 code/MSA_txt_csv.py
```
- aggergate csv files and match input/output pair
```console
python3 code/MAS_input.py
```
Generate processed vectors:
```console
python3 code/data_generation.py
```

train model:

```console
python3 code/training.py
```

run evaluation on the trained models and reproduce figures:

```console
python3 code/evaluations.py
```

