# MAU Simulations
R code to simulate zooarchaeological samples from hypothetical skeletal part profiles.

This repository includes five files. **SEAsim.R** and **SEAexamples.R** are the primary scripts. The former script calls **SEAsim_PlotFunctions.R** and **SEAsim_SimFunctions.R** to create functions necessary to carry out the script. **SEAexamples.R** calls **SEAexamples_Functions.R**, which also generates functions necessary to run the script. As such, the necessary function generating scripts must be in the working directory if running either **SEAsim.R** or **SEAexamples.R**.

**SEAsim.R**

This script simulates draws from a skeletal part profile fossil assemblage. These draws are used to calculate simulated MAU-utility relationships, summarized as Spearman's _rho_, Pearsons's _r_, and Poisson regression _beta_ coefficients. A Bayesian version of the regression coefficients are also generated. Teh code summarizes results as both figures and tables that are written to the working directory.

The script uses three common libraries available in CRAN (_ggplot2_, _parallel_, and _patchwork_) and a fourth package--_Rethinking_, which can be installed from https://github.com/rmcelreath/rethinking.

Lines 22-32 include simulation parameters of interest that may be changed by the user. Note, the "nCo" variable, which specifies the number of samples to draw per sample size ("MNE_vals"), is set to a default value of 100,000. This parameter can have significant impacts on script completion time. On a 36-thread machine, the default value resulted in the simulation finishing in 8 days. If running the script on a personal computer with 8 or fewer threads, it is recommended that the "nCo" argument be set at a value below 100.

**SEAexamples.R**

This script fits Bayesian Poisson regression models to skeletal part data from two archaeological sites and from two experimental skeletal part samples. References for each site and assemblage are specified in the script. It summarizes the results in tables and figures that are written to the working directory. The script should complete in under 15 minutes.
