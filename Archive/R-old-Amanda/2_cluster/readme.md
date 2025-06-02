The cluster code folder contains all of the code for creating antigenic maps on UGA's GACRC. This includes the script for the R environment and the shell script for the GACRC. The code is formulated to allow for OpenMP parallelization of the acmap creation and dimension testing. The output is of the dimension selection criteria and RMSE of all the points.

The Bash code (array_sub.sh) is written to perform an array of the scripts and run them in parallel. The 'input.lst' file contains on each line the names of the datafiles to be used by the Rscript. Each array job will pull one file. The Rscript will then use that Bash variable listed in the submission script as an argument.

If changing the number of files to run the array on:
* change the input.lst file to match the documents
* change the array= variable in the submission script to match the number of arrays

Note: Racmacs will need to be available to the GACRC's system. This can either be done by requesting GACRC downloads the module, or download the module to your personal R library directory ("../workdir/Rlibs/")

To run the script on GACRC add files, and submit the following command:

$ sbatch array_sub.sh

The dimensional analysis takes between 12 to 36 hours for completion.