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
Coinfinder is designed to take as input a dataset of pangenomes and their genes. Ideally, genes will clustered into homologous gene clusters using a pangenomic tool such as <a href="https://github.com/gtonkinhill/panaroo">Panaroo</a>, <a href="https://sanger-pathogens.github.io/Roary/">Roary</a>, <a href="https://github.com/SionBayliss/PIRATE">PIRATE</a>, or <a href="https://github.com/rmcolq/pandora">Pandora</a>. Coinfinder should be used to identify coincident gene sets within a given pangenomic dataset. Coinfinder was written to identify coincident genes among strains of prokaryote species (i.e. a species pangenome) but can be extended to other pangenomic datasets.

### Where can I read more about it? ###
Fiona J. Whelan, Martin Rusilowicz, & James O. McInerney. "<b>Coinfinder: detecting significant associations and dissociations in pangenomes</b>." <a href="https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000338?originator=authorOffprint&identity=174020&timestamp=20210225172138&signature=1155c94b880e4257e2c7bf8f543b7d1c#header">doi: https://doi.org/10.1099/mgen.0.000338</a>

------

### Installation: ###

#### Bioconda ####

If you use Conda: `conda install -c defaults -c bioconda -c conda-forge coinfinder`

(If the installation gets stuck on solving the environment, please see issue <a href="https://github.com/fwhelan/coinfinder/issues/36">36</a>.)

#### Native install ####

##### Dependencies: #####

* `Cmake3.6` or greater https://cmake.org/download/ (`brew install cmake` on a Mac)
* `Python>3.6;<3.8` https://www.python.org/downloads/
* `Boost1.66` or greater https://www.boost.org/users/download/ (`brew install boost` on a Mac)<!-- If there are any issues with Boost, on my new Macbook, I found I had to do brew install boost-python via the below post to get cmake to recognize Boost properly. https://stackoverflow.com/questions/26024878/cmake-cannot-find-boost-on-os-x-with-brew-->
* `OpenMP` (`brew install llvm` on a Mac)<!--Mac: brew install libomp; brew install llvm; -->
* `gcc 5` or greater (default on most new-ish machines)
* `R` libraries: `caper, phytools, getopt, igraph, dplyr, cowplot, data.table, ggraph, flock, future`
* Bionconductor `R` library: `ggtree` https://bioconductor.org/packages/release/bioc/html/ggtree.html

##### Installation: #####

```
cmake -DCMAKE_BUILD_TYPE=Release .
cmake --build .
./coinfinder
```

On macOS, the default compiler may be `clang` instead of `g++`. If so, you may need to point the compiler to gcc; for example:
`export CC=/usr/local/bin/gcc-6; CXX=/usr/local/bin/g++-6; MPICXX=/usr/local/bin/mpic++`

------

### Usage: ###

`coinfinder -i <gene information> [-I] -p <phylogeny> -o <output prefix> [--associate|--dissociate]`

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

Note: the `gene_presence_absence.csv` output from Panaroo appears to differ from Roary in that fields are not surrounded by double-quotes. Coinfinder assumes this double-quote format; you could use something like the following to correct for this:

```
sed -e 's/^/"/g' -e 's/$/"/g' -e 's/,/","/g' gene_presence_absence.csv > gene_presence_absence-withquotes.csv
```

The phylogeny should be Newick-formatted with no zero-length branches. We suggest that this phylogeny be constructed using the core gene information (for example, as suggested in the Roary pipeline https://sanger-pathogens.github.io/Roary/).

Lastly, the user must decide between running Coinfinder to find associations (gene pairs present together) or dissociations (gene pairs which are present apart, or avoid each other).

For more information on usage, please see `coinfinder -h`:

```
File input- specify either: 
    -i or --input          The path to the gene_presence_absence.csv output from Roary
                           -or-
                           The path of the Alpha-to-Beta file with (alpha)(TAB)(beta)
    -I or --inputroary     Set if -i is in the gene_presence_absence.csv format from Roary
    -p or --phylogeny      Phylogeny of Betas in Newick format (required)
Max mode (mandatory for coincidence analysis):
    -a or --associate      Overlap; identify groups that tend to associate/co-occur.
    -d or --dissociate     Separation; identify groups that tend to dissociate/avoid.
Significance- specify: 
    -L or --level          Specify the significnace level cutoff (default: 0.05)
Significance correction- specify: 
    -m or --bonferroni     Bonferroni correction multiple correction (recommeneded)
    -n or --nocorrection   No correction, use value as-is
    -c or --fraction       (Connectivity analysis only) Use fraction rather than p-value
Alternative hypothesis- specify: 
    -g or --greater        Greater (recommended)
    -l or --less           Less
    -t or --twotailed      Two-tailed
Miscellaneous:
    -x or --num_cores      The number of cores to use (default: 2)
    -v or --verbose        Verbose output.
    -r or --filter         Permit filtering of saturated and low-abundance data.
    -U or --upfilthreshold Upper filter threshold for high-abundance data filtering (default: 1.0 i.e. any alpha in >=100/% of betas.
    -F or --filthreshold   Threshold for low-abundance data filtering (default: 0.05 i.e. any alpha in <=5% of betas.
    -q or --query          Query a specific gene.
    -T or --test           Runs the test cases and exits.
    -E or --all            Outputs all results, regardless of significance.
Output:
    -o or --output         The prefix of all output files (default: coincident).
```

### Example output: ###
<div align="center">
<p align="center">
    <img src="Figure1.png?raw=true?" alt="example-output" width="900">
</p>
</div>
An example association network in which each gene (node) is connected to another gene with a line (edge) iff they statistically co-occur with each other. Nodes are weighted by lineage-independence in the phylogeny (i.e. the larger the node, the more phylogenetically independent the gene). Nodes are coloured by connected component, or the set of genes with associative relationships with each other. This data can also be shown as a presence/absence heatmap in relation to the phylogeny (note: this heatmap is a subset of all results; in particular, the large wine coloured gene set has been removed for ease of visibility).


### Example usage: ###

The example dataset, including input and expected output files using the associated manuscript can be found <a href="https://github.com/fwhelan/coinfinder-manuscript">here</a>.

### Citation information: ###
@article{mbs:/content/journal/mgen/10.1099/mgen.0.000338,
   author = "Whelan, Fiona Jane and Rusilowicz, Martin and McInerney, James Oscar",
   title = "Coinfinder: detecting significant associations and dissociations in pangenomes",
   year = "2020",
   publisher = "Microbiology Society",
   url = "https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000338",
   doi = "https://doi.org/10.1099/mgen.0.000338",
   keywords = "pangenome",
   keywords = "gene association networks",
   keywords = "gene co-occurrence",
   abstract = "The accessory genes of prokaryote and eukaryote pangenomes accumulate by horizontal gene transfer, differential gene loss, and the effects of selection and drift. We have developed Coinfinder, a software program that assesses whether sets of homologous genes (gene families) in pangenomes associate or dissociate with each other (i.e. are ‘coincident’) more often than would be expected by chance. Coinfinder employs a user-supplied phylogenetic tree in order to assess the lineage-dependence (i.e. the phylogenetic distribution) of each accessory gene, allowing Coinfinder to focus on coincident gene pairs whose joint presence is not simply because they happened to appear in the same clade, but rather that they tend to appear together more often than expected across the phylogeny. Coinfinder is implemented in C++, Python3 and R and is freely available under the GNU license from <a href="https://github.com/fwhelan/coinfinder" target="xrefwindow">https://github.com/fwhelan/coinfinder.</a>
            "
}


### What if I find a bug or have an issue running coinfinder? ###

If you run into any issues with coinfinder, we want to hear about it! Please don't be shy, and log an Issue including as much of the following as possible:  

* The exact command that you used to call coinfinder (helps us identify where in the code the bug might be).  
* A reproducible example of the issue with a small dataset that you can share (helps us identify whether the issue is specific to a particular computer, operating system, and/or dataset).  
