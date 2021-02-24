# Assessing the significance of directed and multivariate dependence measures
Copyright (C) 2020 Oliver Cliff.

This repository provides MATLAB functions for computing and assessing the linear dependence between multiple autocorrelated time series. This includes various linear dependence measures and the hypothesis tests for inferring their significance, all discussed in our paper in [Phys. Rev. Research](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.013145) and [arXiv](https://arxiv.org/abs/2003.03887).

The measures implemented are: **mutual information**, **conditional mutual information**, **Granger causality**, and **conditional Granger causality** (each for *univariate* and *multivariate* linear-Gaussian processes). For completeness we have also included **Pearson correlation** and **partial correlation** for *univariate* processes (with a potentially *multivariate* conditional process).

The code is licensed under the [GNU GPL v3 license](http://www.gnu.org/licenses/gpl-3.0.html) (or later).

# Getting started
1. Ensure you have MATLAB's [Econometrics Toolbox](https://www.mathworks.com/products/econometrics.html) and [Signal Processing Toolbox](https://www.mathworks.com/products/signal.html) downloaded and installed (for the autocorrelation and filtering functions).
2. Clone (or download) the repository.
3. Add the repository to your path (including the `utils` folder), e.g., by using: ``addpath(genpath('/path/to/repository')).``
6. Documentation found in help for each function. The main functions used are `mvmi.m` (**mutual information**), `mvgc.m` (**Granger causality**), and `pcorr.m` (**Pearson/partial correlation**). All three allow adding a conditional process and can optionally output a *p*-value. The *p*-values are either generated from the finite/asymptotic tests (i.e., the **chi-square** and **F-tests**), or the modified tests (**modified F-** and **lambda-tests**) that we derive in [our paper](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.013145).
7. Demos are included in the `demos` subfolder, including all experiments from the paper. Start with `demos/demoForPaper.m` or `demos/hcpCaseStudy.m` to see all the results from the paper. As an example of modifying Student's t-tests for **Pearson correlation**, we have also included a demo on inferring cross-correlations between *univariate* autocorrelated processes in `demos/demoCrossCorr.m`.

# Citation
Please **cite** your use of this code as:

Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine, Joseph T. Lizier, "Assessing the significance of directed and multivariate measures of linear dependence between time series," Phys. Rev. Research 3, 013145 (2021).

# Acknowledgements
This project was in part supported through:

- The Australian Research Council DECRA grant DE160100630.
- A University of Sydney Robinson Fellowship and NHMRC Project Grant 1156536.
- The University of Sydney Research Accelerator (SOAR) Fellowship program.
- The Centre for Translational Data Science at The University of Sydney’s Research Incubator Funding Scheme.
