library(rstudioapi)
this.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(dirname(this.dir ) )

#setwd("../)
###################################Real data analysis; Get most likeli order (Additional file 2 XLSX)###############################################

source("code/run_yeast.R")

source("code/run_pombe2.R")

source("code/run_fly2.R")

source("code/run_zebrafish2.R")

source("code/run_human.R")

rm(list=ls())
###################################Real data analysis; Analysis of intron splicing order###################################

####Compare in order spliced properties between different species (Figure 3)
source("code/analysis/comp_in_order_p_betw_species.R")

####Get the relative order of the first and last intron in each species (Figure 4)
source("code/analysis/first_last_iso_intron_dis.R")

####correlation intron metric with most likely order (Supplementary figure 2 in Additional file 1 DOCX)
source("code/analysis/cor_cal_call.R")

rm(list=ls())


###################################Simulation analysis######################################################################
####Simulation based on pombe transcriptome (Figure 2A), simulation based on human transcriptome (Figure 2B),
#this depends other software (minimap2,samtools,STAR,Java),their source code are attached in code/simulation_pombe/ and code/simulation_human/

####Simulation the read count matrix directly (Figure 2C)
source("code/analysis/test_model.R")
test_test_large_intron_given_order_times(200);

####Compare hill climbing and permutation (Supplementary figure 5 in Additional file 1 DOCX)
source("code/analysis/test_algorithm.R");

#Generate simulated read count matrix with random read count (Supplementary figure 5A)
scater_compare_simulate(200);
#Fill read count matrix with human real dataset (Supplementary figure 5B)
scater_compare_real_mat();



rm(list=ls())

