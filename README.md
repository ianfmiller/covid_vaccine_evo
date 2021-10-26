# Code for: Assessing the risk of vaccine-driven virulence evolution in SARS-CoV-2
## Analyses
### All analyses presented in the main text are run in the "analysis.R" script. Comments in the code file describe the flow of analyses.
### The analysis script can be sourced to replicate all of the analyses, but the runtime is extremely long. We recommend either running individual analyses on a personal computer (expect runtimes of 2-4hrs), or making use of a computing cluster.
### The epidemiological model is implemented in C and contained in the file "epi.model.c".
### Convenience functions used in the analyses are contained in the "functions.R" script. Comments in the code file describe the usage of the functions.
### An example of setting b1 and b2 parameters is given in the "setting b1 and b2 example.R" file.
## Data visualization
### Output from analyses is stored in "/sim.data"
### Figures are generated in the code files with corresponding names.
