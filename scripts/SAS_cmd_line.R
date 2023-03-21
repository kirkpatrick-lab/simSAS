#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#needs your data file, a window
#and nruns, an int >=1
tot_time_start = Sys.time()
source("scripts/beta_f_sampling.R")
source("scripts/likelihood_function_vectorized.R")
source("scripts/SAS_mod.R")
source("scripts/ML_functions.R")
library(abc)
library(tidyverse)
library(data.table)
library(parallel)

data = args[1]
nruns = args[2] #try 10k to start 

set.seed(123)
#beta_prior = rnorm(nruns, mean = 0.004, sd = .004)
#NEW PRIORS! Acc to marks doc 
x_b = runif(nruns, 0, 1)
beta_prior = 10^(-(4*x_b + 1))

#f_prior = runif(nruns, min = 0, max = 1)
set.seed(123)
#f_prior = rbeta(nruns, 0.001, 1)
x_f = runif(nruns, 0, 1)
f_prior = 10^(-3*x_f)

set.seed(123)
fset = runif(nruns, 0, 1)
beta_prior[fset > f_prior]= 0

data <- read.table(data, header = T) %>% group_by(Window) %>%
  mutate(all_counts = male_counts + female_counts)
nM_byw = data %>% summarise(nM = sum(male_counts))
nF_byw = data %>% summarise(nF = sum(female_counts))
hap_freq <- data %>% mutate(freqF = female_counts/sum(female_counts), freqM = male_counts/sum(male_counts), freq_all = all_counts/sum(all_counts))

#ml_output_t <- lapply(winlist, likelihood.out.vec, data = data, windowsize = 5)
ml_output_t <- likelihood.out.vec(data, windowsize = 5, window = data$Window[1])
#from parallelized 
#ml_output_t <- do.call(rbind, ml_output_t)
ml_output_t$Window <- as.numeric(ml_output_t$Window)

#don't need to order, because all windows will be the same 
#ml_output_t[order(ml_output_t$Window),] -> ml_output_t


out = data.frame()
for(k in 1:nruns){
  res = p_and_beta_from_bf(data = data, nM_byw = nM_byw,
                        nF_byw = nF_byw, hap_freq =  hap_freq, beta = beta_prior[k], f_parm = f_prior[k],
                        ml_output_t)
  out = rbind(out, res)
}


write.table(out, paste(args[1], ".sim", sep = ""), quote = F, row.names = F, col.names = T)
write.table(ml_output_t, paste(args[1], ".true", sep = ""), quote = F, row.names = F, col.names = T)
tot_time_end = Sys.time()
write((paste("Total time for", nruns, "was:", tot_time_end - tot_time_start)), "runlog", append = T)


