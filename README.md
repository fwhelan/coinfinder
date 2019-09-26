# README #

<div align="center">
<p align="center">
    <img src="coinfinder.png?raw=true?" alt="coinfinder-logo" width="200">
</p>
<h1>Coinfinder</h1>
<h3>A tool for the identification of coincident (associating and dissociating) genes in pangenomes.</h3>
Written in collaboration with <a href="https://github.com/mjr129">Martin Rusilowicz</a>.
</div>

### What is it? ###
Coinfinder is an algorithm and software tool that detects genes which associate and dissociate with other genes more often than expected by chance in pangenomes. Coinfinder is written primarily in C++ and is a command line tool which generates text, gexf, and pdf outputs for the user.

Coinfinder uses a Bonferroni-corrected Binomial exact test statistic of the expected and observed rates of gene-gene association to evaluate whether a given gene pair is coincident.

### When and why should I use it? ###
Coinfinder is designed to take as input a dataset of pangenomes and their genes. Ideally, genes will clustered into homologous gene clusters using a pangenomic tool such as <a href="https://sanger-pathogens.github.io/Roary/">Roary</a>, <a href="https://github.com/SionBayliss/PIRATE">PIRATE</a>, or <a href="https://github.com/rmcolq/pandora">Pandora</a>. Coinfinder should be used to identify coincident gene sets within a given pangenomic dataset. Coinfinder was written to identify coincident genes among strains of prokaryote species (i.e. a species pangenome) but can be extended to other pangenomic datasets.


### Dependencies: ###

* `Cmake3.6` or greater https://cmake.org/download/
* `Python3.6` or greater https://www.python.org/downloads/
* `Boost1.66` or greater https://www.boost.org/users/download/
<!-- If there are any issues with Boost, on my new Macbook, I found I had to do brew install boost-python via the below post to get cmake to recognize Boost properly. https://stackoverflow.com/questions/26024878/cmake-cannot-find-boost-on-os-x-with-brew-->
* `OpenMP` 
* `R` libraries: `caper, phytools, getopt, igraph, dplyr, cowplot, data.table, ggraph`
* Bionconductor `R` library: `ggtree` https://bioconductor.org/packages/release/bioc/html/ggtree.html

### Quick installation instructions: ###

```
cmake -DCMAKE_BUILD_TYPE=Release .
cmake --build .
./coinfinder
```

### Usage: ###

`coinfinder [-d|-D] <gene information> -p <phylogeny> -o <output prefix> [--associate|--dissociate]`

Coinfinder requires gene information and a phylogeny as input. The gene information can be provided in one of two formats: (a) as the `gene_presence_absence.csv` output from <a href="https://sanger-pathogens.github.io/Roary/">Roary</a>; (b) as a tab-delimited list of genes present in each strain. An example of a tab-delimited list of genes:

```
gene_1	genome_1
gene_1	genome_2
gene_1	genome_3
gene_2	genome_2
gene_2	genome_3
gene_3	genome_1
gene_3	genome_2
```

The phylogeny should be Newick-formatted with no zero-length branches. We suggest that this phylogeny be constructed using the core gene information (for example, as suggested in the Roary pipeline https://sanger-pathogens.github.io/Roary/).

Lastly, the user must decide between running Coinfinder to find associations (gene pairs present together) or dissociations (gene pairs which are present apart, or avoid each other).

### Example usage: ###

Coming soon...

<!--### Citation information (pre-print): ###
@article {}-->


### What if I find a bug or have an issue running coinfinder? ###

If you run into any issues with coinfinder, we want to hear about it! Please don't be shy, and log an Issue including as much of the following as possible:  

* The exact command that you used to call coinfinder (helps us identify where in the code the bug might be).  
* A reproducible example of the issue with a small dataset that you can share (helps us identify whether the issue is specific to a particular computer, operating system, and/or dataset).  
