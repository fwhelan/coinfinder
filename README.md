# README #

<div align="center">
<p align="center">
    <img src="coinfinder.png?raw=true?" alt="coinfinder-logo" width="200">
</p>
Coinfinder
A tool for the identification of coincident (co-occurring and avoiding) genes in pangenomes.
Written in collaboration with <a href="https://github.com/mjr129">Martin Rusilowicz</a>.
</div>

### Quick installation instructions: ###

```
cmake -DCMAKE_BUILD_TYPE=Release .  
cmake --build .  
./coinfinder  
```

### Dependencies: ###

* `Cmake3.6` or greater
* `Python3.6` or greater
* `Boost1.66` or greater 
* `OpenMP`
* `R` libraries: `caper, phytools, getopt, igraph, dplyr, cowplot, data.table, ggtree` (from Bioconductor), `ggraph`

### Example usage: ###

Coming soon...  


### What if I find a bug or have an issue running coinfinder? ###

If you run into any issues with coinfinder, we want to hear about it! Please log an Issue using the toolbar to the left. In the Issue report please include as many of the following as possible:  

* The exact command that you used to call coinfinder (helps us identify the parts of the code that the bug might reside in).  
* A reproducible example of the issue with a small dataset that you can share (helps us identify whether the issue is specific to a particular computer and/or operating system).  
