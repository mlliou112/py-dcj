# Py-DCJ

## What is Py-DCJ?

Py-DCJ is a python implementation of a distance algorithm based on the double-cut and join (DCJ) model of genome rearrangements. DCJ is an abstract operation on sets of linear and circular chromosomes defined by cutting the genome in 2 locations and rejoining the ends. General information about the model can be found on the [wikipedia page](https://en.wikipedia.org/wiki/Double_Cut_and_Join_Model).


## Features

Chromosome representation under the DCJ model is a list of unique markers representing "genes" or a block of dna. Linear chromosomes are denoted with a telomere `"."` at each end of the chromosomes, while circular chromosomes do not have telomeres. Genomes are represented as lists of chromosomes. 

For example `[".abcD.", "eF"]` is a genome with 1 linear and 1 circular chromosome. Each unique marker is a single letter. Capital letters indicate that the gene is in reverse orientation. 

Py-DCJ can calculate distances based on common markers between the two genomes (default method: dcj) , or including insertion-deletion (indel) information. For more details about the distance algorithm, see [this paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.911.7351&rep=rep1&type=pdf).

## Example
```
>>> from dcj.main import calculate_distance 
>>> genome_a = [".abc."]
>>> genome_b = [".aCBd."]
>>> calculate_distance(genome_a, genome_b, method="dcj")
... 1.0
>>> caculate_distance(genome_a, genome_b, method="indel")
... 2.0
```

## Installation

The recommended way to install Py-DCJ is via git clone.

