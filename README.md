# Developmental Trajectories

Characterizing the early-life gut microbiome using 16S data from the COPSAC2010 Cohort

## Scripts

This repository contains the R scripts nessecary to recreate plots in my masters thesis submited 21st of May 2024.

- packages.R: Installs the packages necessary to run the remaining scripts
- alpha.R: Characterizes alpha diversity and depth biases using rarefaction
- beta.R: Characterizes beta diversity using rarefaction
- tcgsa.R: A readable version of the orginal time-course gene set analysis script
- assign.R: An alternative assignment algorthim
- bubblegum.R: Bubblegum plot to chracterize the impact of environmental factors

## Data

Unfortunately, I cannot share the files here as they are confidential. Contact me or COPSAC to get them

- dada2_gtdb_2021_11_18.RData: Output phyloseq from the dada2 pipeline. ASV matrix is counts
- Feaces_16S_Traj_DMM.RData: Output from the original TcGSA script
- exposure5.RData: Metadata on 5 key enviromental factors
- bmi_6y.RData: Metadata on a lot of environmental factors
- dist_rarefaction_100_2020.RData: Bray Curtis distance matrix, output from beta.R
- phy_assigned.RData: Abundance is mean relative abundance transformed, output from assign.R

## Contact

Anton Kjellberg <br/>
antonkjellberg@hotmail.com <br/>
+45 42 22 64 13 <br/>
Collegium Domus Regiæ <br/>
St. Kannikestræde 2 <br/>
1169 Kbh. K <br/>



