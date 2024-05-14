install.packages('tidyverse') # Data processing
install.packages('patchwork') # Combine plots
install.packages('rstatix') # Pipe-friendly statistic tests
install.packages('ggpubr') # # Statistics and tables in plots
install.packages('ggridges') # Ridgeline plots
install.packages('ggside') # Plots on the sides of plots
install.packages('gghalves') # Half-half plots
install.packages('TcGSA') # Time-course gene set analysis
install.packages('BiocManager') # Bioinformatics package installation
BiocManager::install('phyloseq') # Integrated metagenomics objects
BiocManager::install('multtest') # Required for ranacapa and TcGSA
install.packages('devtools') # Install packages from Github
devtools::install_github("mikemc/speedyseq") # Faster phyloseq functions
devtools::install_github("gauravsk/ranacapa") # Fast rarefaction curves
devtools::install_github("kylebittinger/usedist") # Distance matrix utilities