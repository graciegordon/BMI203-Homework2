# Homework 2: Clustering

[![Build
Status](https://travis-ci.org/graciegordon/BMI203-Homework2.svg?branch=master)](https://travis-ci.org/graciegordon/BMI203-Homework2)

Skeleton for clustering project.

## assignment

1. Implement a similarity metric
2. Implement a clustering method based on a partitioning algorithm
3. Implement a clustering method based on a hierarchical algorithm
4. Answer the questions given in the homework assignment


## structure

The main file that you will need to modify is `cluster.py` and the corresponding `test_cluster.py`. `utils.py` contains helpful classes that you can use to represent Active Sites. `io.py` contains some reading and writing files for interacting with PDB files and writing out cluster info.

```
.
├── README.md
├── data
│   ...
├── functions
│   ├── __init__.py
│   ├── __main__.py
│   ├── cluster.py
│   ├── io.py
│   └── utils.py
└── test
    ├── test_cluster.py
    └── test_io.py
```

## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `functions/__main__.py`) can be run as
follows

```
python -m functions -P data test.txt
```
The -P flag is for the partition method, the -H flag is for heirarchical,
and the -A flag runs both methods, then compares the two clusters using the
Rand Index and measures the Quality of the clusters using the Sum of Distances
Method

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.
