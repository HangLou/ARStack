ARStack
========================================

## Environment Setup
The code has been tested successfully using Python 3.8; thus we suggest using this version or a later version of Python. A typical process for installing the package dependencies involves creating a new Python virtual environment.

To install the required packages, run the following:
```console
pip install .
```
## running codes
Generate processed vectors:
```console
python3 code/data_generation.py
```

train model:

```console
python3 code/training.py
```

run evaluation on the trained models:

```console
python3 code/evaluations.py
```

Reproduce plots:
```console
python3 code/plots.py
```

## TODO:
1. Access to the raw data of all 3 dataset.
2. Code and data to produce the MSA features csv.
3. Add functions in code/plots.py for the plots in the paper. 