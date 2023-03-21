ssd_from_beta_f <- function(data, nM_byw, nF_byw, hap_freq, beta, f_parm, ml_output_t){
    #data is the data frame we read in at the beginning describing the observed data 
    #beta and F are derived from uniform distributions, but modified as per email from Mark
    
    uniq_winds = unique(data$Window)
    n = 5
    k = 3
    num = 1
    
    allout = data.frame()
    for(window in uniq_winds){
      winout = data %>% filter(Window == window)
      samp_m = nM_byw %>% filter(Window == window) %>% select(nM) %>% as.numeric()
      samp_f = nF_byw %>% filter(Window == window) %>% select(nF) %>% as.numeric()
      #samp = sum(winout$all_counts)
      p = winout$Center_MAF[1]
      f = hap_freq %>% filter(Window == window) %>% ungroup() %>% select(freq_all) %>% unlist() %>% as.vector()
      sim1 = simSAS(p, f, n, k, beta, samp_m, samp_f, num)
      winout$female_counts = sim1$lnfemales[,1]
      winout$male_counts = sim1$lnmales[,1]
      winout$all_counts = winout$female_counts + winout$male_counts
      allout = rbind(allout, winout)
    }

    #now run jared's pipeline on fake1 and window1
    #inefficient to run this over every loop iteration - do it in a separate
    #function 
    #the likelihood for the original data can be calculated outside of the function - 

    winlist <- as.character(unique(allout$Window))
    #artificially inflate by 1
    allout$female_counts = allout$female_counts + 1
    allout$male_counts = allout$male_counts + 1
    ml_output_f <- mclapply(winlist, likelihood.out.vec, data = setDT(allout), windowsize = 5) #can't handle zeroes
    ml_output_f <- do.call(rbind, ml_output_f)
    ml_output_f$Window <- as.numeric(ml_output_f$Window)
    ml_output_f[order(ml_output_f$Window),] -> ml_output_f
    
    ssd = sum((sort(ml_output_t$p_value) - sort(ml_output_f$p_value))^2)
    
    #ml_output_f #why do my nF0 and nM0s not match the amount in the dataset - there must be a posthoc calculation here 
    
    #now let's get the sum of squared difference between them
    #percentile_vals =  set a percentile with the real data 
    return(cbind(beta, f_parm, ssd))
}

p_and_beta_from_bf <- function(data, nM_byw, nF_byw, hap_freq, beta, f_parm, ml_output_t){
  #data is the data frame we read in at the beginning describing the observed data 
  #beta and F are derived from uniform distributions, but modified as per email from Mark
  
  uniq_winds = unique(data$Window)
  n = 5
  k = 3
  num = 1
  
  allout = data.frame()
  #for(window in uniq_winds){
  #non parallelized for windows
    window = uniq_winds
    winout = data %>% filter(Window == window) #reduce this later, in theory should only have one window per input file 
    samp_m = nM_byw %>% filter(Window == window) %>% select(nM) %>% as.numeric()
    samp_f = nF_byw %>% filter(Window == window) %>% select(nF) %>% as.numeric()
    #samp = sum(winout$all_counts)
    p = winout$Center_MAF[1]
    f = hap_freq %>% filter(Window == window) %>% ungroup() %>% select(freq_all) %>% unlist() %>% as.vector()
    sim1 = simSAS(p, f, n, k, beta, samp_m, samp_f, num)
    winout$female_counts = sim1$lnfemales[,1]
    winout$male_counts = sim1$lnmales[,1]
    winout$all_counts = winout$female_counts + winout$male_counts
    #allout = rbind(allout, winout)
  #}
  
  #now run jared's pipeline on fake1 and window1
  #inefficient to run this over every loop iteration - do it in a separate
  #function 
  #the likelihood for the original data can be calculated outside of the function - 
  allout = winout #sloppy, fix later 
  winlist <- as.character(unique(allout$Window))
  #artificially inflate by 1
  allout$female_counts = allout$female_counts + 1
  allout$male_counts = allout$male_counts + 1
  #ml_output_f <- mclapply(winlist, likelihood.out.vec, data = setDT(allout), windowsize = 5) #can't handle zeroes
  #do not in parallel
  ml_output_f <- likelihood.out.vec(data = allout, windowsize = 5, window = data$Window[1])
  #ml_output_f <- do.call(rbind, ml_output_f)
  ml_output_f$Window <- as.numeric(ml_output_f$Window)
  #ml_output_f[order(ml_output_f$Window),] -> ml_output_f
  
  #ssd = sum((sort(ml_output_t$p_value) - sort(ml_output_f$p_value))^2)
  
  #ml_output_f #why do my nF0 and nM0s not match the amount in the dataset - there must be a posthoc calculation here 
  
  #now let's get the sum of squared difference between them
  #percentile_vals =  set a percentile with the real data 
  output = cbind(beta, f_parm, ml_output_f)
  colnames(output) = c("beta_sim", "f_sim", colnames(ml_output_f))
  return(output)
}

