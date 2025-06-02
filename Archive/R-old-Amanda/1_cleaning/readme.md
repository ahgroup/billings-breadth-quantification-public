Packages needed for running the code:

Biostrings, devtools, here, Rcpi, Racmacs, seqinr, tidyverse, ggrepel, ggpubr, cowplot, latex2exp

To install Biostrings package: https://bioconductor.org/packages/release/bioc/html/Biostrings.html
To install Rcpi package: https://bioconductor.org/packages/release/bioc/html/Rcpi.html
To install Racmacs package: https://github.com/acorg/Racmacs
* Be sure to read through the introduction to confirm that it is configured correctly. Follow the site for troubleshooting details for download, compiling, and function.


# Files

1_cleaning_code.R

* This file pulls in the raw datasets and creates processed data. It also creates the five subsetted datasets.

2_pepitope_calculator.R

* This file calculates the P-epitope values for the pairwise distances. There are save toggles. They are set to save files and figures by default. There is a save extra toggle that is set to FALSE. By default the difference matrices are not saved. This can be set to TRUE to save the distance matrices. 

3_acmap_creation.R

* Creates 10 maps for the subset analysis, 5 for H1N1, 5 for H3N2 (SD only, post-vaccination)
* Creates 12 maps per subtype (H1N1, H3N2), for: SD pretiter 1, 2, and 3 dimensions (3 maps), SD posttiter for 1, 2, 3 dimensions (3 maps), all individuals pretiter 1, 2, 3 dimensions (3 maps), and all individuals posttiter 1, 2, 3 dimensions (3 maps)
* Creates the 8 input files for the dimension analysis on the high computing cluster (H1/H3 by Pre/Post by All/SD = 8 total). If user doesn't have access to the a high performance computing cluster, the output files are saved as well so that the step can be skipped if necessary. 
* Creates the 2 files for the distances calculated from the maps created from each individual season rather than all of the sera from all of the seasons together on one map

4_acmap_distance.R

4_acmap_subsets.R

5_fonville_datasets.Rmd