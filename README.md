# Code for: Assessing the risk of vaccine-driven virulence evolution in SARS-CoV-2
## Analyses
### All analyses presented in the main text are run in the "analysis.R" script. Comments in the code file describe the flow of analyses.
### The analysis script can be sourced to replicate all of the analyses, but the runtime is extremely long. 
### Analyses can alternatively be run on a SLURM based computer cluster by executing the 'batch.script.q' script
#### Each analysis is split into it's own directory
#### Directories 1-12 correspond to the main text analyses.
#### Directories 13-21 correspond to stronger effects of natural immunity than assumed in the main text (presented in Fig. S4)
#### Directories 22-30 correspond to weaker effects of natural immunity than assumed in the main text (presented in Fig. S5)
#### Directories 22-30 correspond to smaller contributions of the LRT to transmission than assumed in the main text (presented in Fig. S3)

We recommend either running individual analyses on a personal computer (expect runtimes of 2-4hrs), or making use of a computing cluster.
### The epidemiological model is implemented in C and contained in the file "epi.model.c".
### Convenience functions used in the analyses are contained in the "functions.R" script. Comments in the code file describe the usage of the functions.
### An example of setting b1 and b2 parameters is given in the "setting b1 and b2 example.R" file.
## Data visualization
### Output from analyses is stored in "/sim.data"
### Figures are generated in the code files with corresponding names.
