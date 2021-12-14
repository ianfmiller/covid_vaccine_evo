# Code for: Assessing the risk of vaccine-driven virulence evolution in SARS-CoV-2
## Analyses
All analyses presented in the main text are run in the "analysis.R" script. Comments in the code file describe the flow of analyses.
As the runtime of the "analysis.R" script is extremely long, analyses can be run in parallel on a SLURM based computer cluster by executing the "batch.script.q" script in the "cluster" file. This script runs individual analyses located in the directories "dir.1" through "dir.39".
Directories "dir.1"-"dir.12" correspond to the main text analyses.
Directories "dir.13"-"dir.21" correspond to stronger effects of natural immunity than assumed in the main text (presented in Fig. S5).
Directories "dir.22"-"dir.30" correspond to weaker effects of natural immunity than assumed in the main text (presented in Fig. S4).
Directories "dir.31"-"dir.39" correspond to smaller contributions of the LRT to transmission than assumed in the main text (presented in Fig. S3).
The epidemiological model is implemented in C and contained in the file "epi.model.c".
Convenience functions used in the analyses are contained in the "functions.R" script. Comments in the code file describe the usage of the functions.
An example of setting b1 and b2 parameters is given in the "setting b1 and b2 example.R" file.
## Data visualization
Output from analyses is stored in "/sim.data".
Figures are generated in the code files with corresponding names.
