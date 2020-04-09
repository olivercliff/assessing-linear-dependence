# Exact Inference of Linear Dependence
Copyright (C) 2020 Oliver Cliff.

This repository provides MATLAB functions for the exact inference of linear dependence between multiple autocorrelated time series. This includes various linear-dependence measures and the hypothesis tests for inferring their significance, all discussed in the paper found here:

[https://arxiv.org/abs/2003.03887](https://arxiv.org/abs/2003.03887)

The measures implemented are: **mutual information**, **conditional mutual information**, **Granger causality**, and **conditional Granger causality** (each for *univariate* and *multivariate* linear-Gaussian processes).

The code is licensed under the [GNU GPL v3 license](http://www.gnu.org/licenses/gpl-3.0.html) (or later).

# Getting started
1. Ensure you have MATLAB's [Econometrics Toolbox](https://www.mathworks.com/products/econometrics.html) and [Signal Processing Toolbox](https://www.mathworks.com/products/signal.html) downloaded and installed (for the autocorrelation and filtering functions).
2. Clone (or download) the repository.
2. Documentation found in help for each function. The main functions used are *mvmi.m* (*multivariate mutual information*) and *mvgc.m* (*multivariate Granger causality*). Both allow adding a conditional process and use the *significance.m* function to generate *p*-values (for both the **likelihood-ratio tests** and **exact tests**).
3. Demos are included in the *demos* subfolder, including all experiments from the paper. Start with *demos/numerical_evaluation.m* or *demos/hcp_case_study.m* to see typical use cases of this package.

# Citation
Please **cite** your use of this code as:

Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine, Joseph T. Lizier, "Exact Inference of Linear Dependence for Multiple Autocorrelated Time Series," arXiv preprint arXiv:2003.03887 (2020).

# Acknowledgements
This project was in part supported through:

- The Australian Research Council DECRA grant DE160100630.
- A University of Sydney Robinson Fellowship and NHMRC Project Grant 1156536.
- The University of Sydney Research Accelerator (SOAR) Fellowship program.
- The Centre for Translational Data Science at The University of Sydney’s Research Incubator Funding Scheme.
