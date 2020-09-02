# PEidos

**P**ython writes **Eidos** is a package that allows you to automate the writing of Eidos/SLiM scripts under very constrained specifications. 

It is written to take a nwk tree in as input and output a SLiM script that replicates the tree's topology.

This Newick formatted tree represents a tree-like population demography where generation time is specified in branch lengths. Leaves should be named under SLiM's requirements (p[0-9]\*), internal nodes should be named arbitrarily. 

Pulse-like and continuous admixture events can be specified (these can happen once populations have achieved tree topology in the simulation).

So far, only neutral scenarios are functinal. I plan on adding advantageous mutations and restart conditional on adaptive mutation loss.


## Eidos is already a scripting language. Y U DO DIS?

My motivation to write this package is that SLiM scripts are not so flexible (at least not as Python). Also, my interest lies in these hierarchical scenarios where typing these demographies in a SLiM script can be error-inducing and time-consuming.

