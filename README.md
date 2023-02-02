# meltR.A.paper

Raw data and analysis code for:

"MeltR Software Provides Facile Determination of Biopolymer Thermodynamics from Melting Data"

Jacob P. Sieg1,2,*, Sebastian J. Arteaga3, Brent M. Znosko3, Philip C. Bevilacqua1,2,4,*

1Department of Chemistry, Pennsylvania State University, University Park, PA 16802.
2Center for RNA Molecular Biology, Pennsylvania State University, University Park, PA 16802.
3Department of Chemistry, Saint Louis University, Saint Louis, MO 63103.
4Department of Biochemistry and Molecular Biology, Pennsylvania State University, University Park, PA 16802.

*Correspondence should be directed to jus841@psu.edu and pcb5@psu.edu 

### Cite 

### Contents

  1. Drafts: Directory containing relavant versions of the manuscript
  2. Figures: Directory containing scripts, intermediate versions, and final versions for each figure
  3. R: Directory containing one script, used to document the raw data
  4. SI_files: Directory containing supplemental files
  5. Tables: Directory containing scripts, intermediate versions, and final versions for each table
  6. data-raw: Directory containing all of the raw data, and a script for compiling the raw data into a single data set
  7. data: Directory containing raw data in *.rda format 
  8. man: Documentation for raw data in *.Rd format
  
 ### Instructions for use with Rstudio
 
 1. Clone (Download the repository). Either click "Code" in the uper right hand corner then "Download zip" and unzip, or type the following in a bash terminal:
 
 ```{r}
 git clone https://github.com/JPSieg/meltR.a.paper
 ```
 
 2. Open the meltR.A.paper.Rproj Rproject file with Rstudio.
 
 3. To expose the data to your memory, type the following in your R console:
 
 ```{r}
 devtools::load_all()
 ```
 
 4. One can run scripts by opening them in R studio and running them line by line, as long as you don't change the project directory.
 
 5. One can also remake figures in a Bash terminal (For example the one you might use inside R studio) from the main directory by running.
 
 ```{r}
 #Using Figure 3 as an example
 Rscript Figures/Figure_3_check_meltR.A/Meltwin_comparison_script.R
 ```

Note, paths are relative so scripts that require reading a data file will not work outside of this dirrectory.
