# LRA

[![Join the chat at https://gitter.im/keyavi/LRA](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/keyavi/LRA?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
Latent RNA-seq Analysis

chat about the project here
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/keyavi/LRA?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

This project explores some methods for clustering of RNA sequence reads without a reference transcriptome.

## Motivation

The goal of this clustering is to reduce the computational complexity of building the de Bruijn Graph. Instead of construction the graph on all sequence reads, the graph would be built separately for each cluster

## Using/Browsing the Code

The source/ directory contains all of the code.

### Locality Sensitive Hashing (LSH)

    lsh.py
    lsh_functions.py

Example of running LSH:

    ./lsh.py -s 10 -k 21 r1.fastq r2.fastq

This runs with k-mer lengths of 21 and a hash size on 10.


### Three Stage Clustering (TSC)

    ThreeStageClustering.py
    Clustering.py
    Cluster.py

Example of running TSC:

    ./ThreeStageClustering.py r1.fastq r2.fastq
