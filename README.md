## bayespetr : Node-wise Pseudo-marginal Algorithm with vSMC Rcpp functionality


### Summary

This package provides R with access to the [NWPM algorithm](https://arxiv.org/abs/2109.08573)
with Rcpp bindings for [vSMC](https://github.com/zhouyan/vSMC) by Zhou(Journal of Statistical Software, 2015, v62, i9)

A toy example, from the above paper, is provided.
This package requires the vSMC C++ library to perform PET image analysis using the NWPM algorithm.
To reproduce tests, replace the files in director , with the files provided in /vSMC_R_rng/ given in bayespetr.

Further integration and extensions are planned.

### PET Image Analysis
A look-up table of $C_T$ will also need to be generated.

### Help
For support and discussion please contact me.

### Author

Denish Thesingarajah
