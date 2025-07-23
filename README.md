On the inclusion of non-concurrent controls in platform trials with an
interim analysis
================

This is an accompanying repository for the paper *“On the inclusion of
non-concurrent controls in platform trials with an interim analysis’’*
by Pavla Krotka, Martin Posch, and Marta Bofill Roig. It contains the
code to reproduce all simulations and figures, as well as the case study
presented in the paper.
<!-- "[On the inclusion of non-concurrent controls in platform trials with an interim analysis](https://arxiv.org)". -->

The repository is structured as follows:

- Folder **functions**:

This folder contains all functions necessary to run the simulation
study:

    - `runtrial_IA_MAE.R` - function for running a platform trial with 2 treatment arms and an interim analysis of arm 1 at the time of adding the arm 2, where arm 2 is evaluated using the mean adjusted estimator
    - `get_mae.R` - function for computing the mean adjusted estimator
    - `bootstrap_functions.R` - function for performing stratified bootstrap with interim analysis to compute the variance of the mean adjusted estimator
    - `sim_function_MAE_par.R` - wrapper function for performing simulation studies for a given set of scenarios

- Folder **paper_sim**:

  - *simscript_paper_final.R*: This script contains the code to
    reproduce all simulations included in the paper and supplementary
    material. The results for the individual scenarios are then saved in
    the subfolder *results_final*.
  - *figures_paper_final.Rmd*: This file contains the code to create all
    figures presented in the paper and the supplementary material. The
    figures are saved in the subfolder *figures* in a .pdf and .tiff
    formats.

## Working directories

The required working directory for each code file is the folder where
this file is located, i.e., for running the simulations, the working
directory should be set to the folder **paper_sim**.

------------------------------------------------------------------------

**Funding**

This publication is supported by the predoctoral program AGAUR-FI ajuts
(2024 FI-1 00401) Joan Oró, which is backed by the Secretariat of
Universities and Research of the Department of Research and Universities
of the Generalitat of Catalonia, as well as the European Social Plus
Fund.

This work was supported by Grant PID2023-148033OB-C21 funded by
MICIU/AEI/10.13039/501100011033 and by FEDER/UE; and the Departament
d’Empresa i Coneixement de la Generalitat de Catalunya (Spain) under
Grant 2021 SGR 01421 (GRBIO).

Marta Bofill Roig is a Serra Húnter Fellow.
