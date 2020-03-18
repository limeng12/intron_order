library(rstudioapi)
this.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) 
#this.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) 

setwd(dirname(getwd() ) )

source("code/run_human.R")

source("code/run_zebrafish2.R")

source("code/run_fly2.R")

source("code/run_pombe2.R")
