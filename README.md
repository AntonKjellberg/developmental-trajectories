# Developmental Trajectories

Characterizing the early-life gut microbiome using 16S data from the COPSAC2010 Cohort

## Bin

This repository contains the R scripts nessecary to recreate plots in my masters thesis submited 21st of May 2024.

- packages.R: Installs the packages necessary to run the scripts
- alpha.R: Alpha diversity and depth biases using rarefaction
- beta.R: Beta diversity using rarefaction
- tcgsa.R: A readable version of the orginal time-course gene set analysis script
- assign.R: An alternative assignment algorthim
- exposures.R: Chracterizes the impact of environmental factors
- scoring.R: Score children based on their deviation and fluctuation compared to mean trajectories

## Data

The data is sensitive information and is not provided here, but it can be acquired through COPSAC.

- dada2_gtdb_2021_11_18.RData: Output phyloseq from the dada2 pipeline. ASV matrix is counts
- exposures.RData: Sample metadata on key enviromental factors
- dist_rarefaction_100_2020.RData: Bray Curtis distance matrix, output from beta.R
- phy_assigned.RData: Cleaned and assigned phyloseq. ASV matrix is counts. Output from assign.R

## Contact

Anton Kjellberg <br/>
antonkjellberg@hotmail.com <br/>
+45 42 22 64 13 <br/>
Collegium Domus Regiæ <br/>
St. Kannikestræde 2 <br/>
1169 Kbh. K <br/>



