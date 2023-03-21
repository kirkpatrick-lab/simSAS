# Functions used in ML pipeline, updated 2/20/22

#Get overlapping windows
getWindows <- function(vec, winsize, overlap) {
  starts = seq(1, length(vec), by=winsize-overlap)
  ends   = starts + winsize - 1
  ends[ends > length(vec)] = length(vec)
  windows <- lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
  names(windows) <- starts
  return(windows)
}

#print windows for haplotype block
window_print <- function(data){ 
  haplotype_counts_boot <- bootstrap(data)
  winds <- NULL
  
  #Likelihood calculations
  winlist <- names(wins)
  wind_out <- lapply(winlist, likelihood.out.boot, data = haplotype_counts_boot, windowsize = win)
  wind_out <- do.call(rbind, wind_out)
  wind_out$Locus <- gsub('\\s+', '', wind_out$Locus)
  wind_out$Window <- gsub('\\s+', '', wind_out$Window)
  wind_out$Window <- as.numeric(wind_out$Window)
  wind_out$Locus <- as.numeric(wind_out$Locus)
  wind_out[order(wind_out$Window, wind_out$Locus),] -> wind_out
  wind_out %>%
    select(Window,Locus) -> wind_out
  winds <- cbind(winds, wind_out$Window)
  winds <- cbind(winds, wind_out$Locus)
  winds <- as.matrix(t(winds))
  cat(winds[1,], "\n")
  cat(winds[2,], "\n")
}

#Get haplotypes
gethaps <- function(wnd, input, s) {
  if (win == 1){
    results = input[input$sex == s, wnd]
  }
  else {
    results = do.call(paste0, input[input$sex == s, wnd])
  }
  return(results)
}

#Counting function
counting <- function(x) {
  u <- unique(x);
  data.frame(
    haplotype=u,
    count=sapply(u, function(v) { length(which(x==v)) } )
  )
}

#Count haplotypes (binary to base 10)
haplocounts <- function(input, sex) {
  convert = lapply(input, strtoi, base="2")
  counts = lapply(convert, counting)
  return(counts)
}


# Determine variants in haplotypes
allele <- function(i,j) {
  x = i/(2^(j-1))
  if (floor(x) %% 2 != 0) {
    return(1)
  } else {
    return(0)
  }
}

#Total haplotypes
total_haps <- function(data,window,hap){
  data %>%
    filter(Window == window, grepl(hap, haplotype)) %>%
    select(counts) %>% sum() -> a
  data %>%
    filter(Window == window) %>%
    select(counts) %>% sum() -> b
  return(a/b)
}

#Sex averaged frequency
sex_avg_frq <- function(data,window,hap){ 
  data %>% filter(Window == window) -> win_hap
  win_hap %>% filter(grepl("f", haplotype)) %>% select(counts) %>% sum() -> sumF
  win_hap %>% filter(grepl("m", haplotype)) %>% select(counts) %>% sum() -> sumM
  win_hap %>% filter(grepl(paste("\\b",hap,"\\b", sep=""), haplotype)) %>% filter(grepl("f", haplotype)) %>% select(counts) %>% unlist() %>% sum(na.rm = TRUE) -> hapF
  win_hap %>% filter(grepl(paste("\\b",hap,"\\b", sep=""), haplotype)) %>% filter(grepl("m", haplotype)) %>% select(counts) %>% unlist() %>% sum(na.rm = TRUE) -> hapM
  a <- hapF/sumF
  b <- hapM/sumM   
  frq <- (a+b)/2
  return(frq)
}

#Haplotype numbers by sex
haplonum <- function(data,hap,window,sex){ 
  data %>% filter(Window == window) -> win_hap
  win_hap %>% filter(grepl(sex, haplotype) & grepl(paste("\\b",hap,"\\b", sep=""), haplotype)) %>%
    select(counts) %>% unlist() %>% sum(na.rm = TRUE)
}

#Check for polymorphism within the haplotypes (eg,for bootstrap sampling)
polycheck_haps <- function(l) {   
  m <- nchar(l[[1L]][1L])
  n <- length(l)
  f0 <- function(x) {
    matrix(unlist(strsplit(x, ""), FALSE, FALSE), m)
  }
  X <- do.call(rbind, lapply(l, f0))
  matrix(matrixStats::rowAnys(X != X[, 1L]), n, byrow = TRUE)
}

# Likelihood function
lik.function <- function(B, params, freqs){ # n1=nF0, n2= nM0, n3= nF1, n4 = nM1, c = c
  n1 <- params[1]
  n2 <- params[2]
  n3 <- params[3]
  n4 <- params[4]
  c <- params[5] #what is c? a constant from the allele freqs?
  p <- freqs[1]     # freqs = p, q
  q <- freqs[2]
  logl <- (n1 * log(1-p*B)) + (n2 * log(1+p*B)) +
    (n3 * log(1+q*B)) + (n4 * log(1-q*B)) + c
  return(logl)
}


#Optimization function site-by-site
likelihood.out <- function(data,windowsize,window){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- data %>% separate(haplotype, into = c("gen", "haplo"), "\\.") %>% select(haplo) %>% mutate(haplo = as.numeric(haplo)) %>% unique()
  haplist = as.vector(haplist$haplo)
  #assign values
  for (s in 1:windowsize){    #site by site
    
    #variables
    frq <- NULL
    numM <- NULL
    numF <-NULL
    h <- NULL
    f <- NULL
    nM_0 <- NULL
    nF_0 <- NULL
    
    for (i in haplist){   
      frq[i+1] <- c(sex_avg_frq(data, window,i))
      numM[i+1] <- c(haplonum(data, i, window, "m"))
      numF[i+1] <- c(haplonum(data, i, window, "f"))
      
      if (frq[i+1] > 0){
        h[i+1] <- c(((numM[i+1]+numF[i+1])*log(frq[i+1])))
      }
      k = ((windowsize+1) - s)
      if (allele(i,k) == 0){  
        f[i+1] <- frq[i+1] 
        nM_0[i+1] <- numM[i+1] 
        nF_0[i+1] <- numF[i+1]
      }
      
    }
    
    #Calculate variables
    nM0 = sum(nM_0, na.rm=TRUE)
    nF0 = sum(nF_0, na.rm=TRUE)
    nM1 = (sum(numM, na.rm=TRUE) - sum(nM0, na.rm=TRUE))
    nF1 = (sum(numF, na.rm=TRUE) - sum(nF0, na.rm=TRUE))
    c = sum(h, na.rm=TRUE)
    loc = colnames(hap_test[wins[[window]][s]]) #what is this?
    hap_N = nrow(samples)*2
    q = ((nM0/(nM0+nM1)) + (nF0/(nF0+nF1)))/2
    p = 1-q
    
    #Optimize the likelihood
    hap_counts <- c(nF0,nM0,nF1,nM1,c)
    allele_freqs <-c(p,q)
    
    #Dynamic bounds for LL search
    b1 <- 0.1
    bL <- -(b1)   #lower bound
    bU <- b1    #upper bound
    e <- 10^(-4)
    
    ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
    
    while ("NaN" %in% ml$objective) {
      b1 <- b1*0.25
      bL <- -(b1)   #lower bound (adjusted)
      bU <- b1    #upper bound (adjusted)
      ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
    }
    
    cL <- abs((ml$maximum) - bL) > e
    cU <- abs((ml$maximum) - bU) > e
    
    while ( cL == "FALSE" | cU == "FALSE"){
      if (cL == "FALSE"){
        bL <- bL - b1
        bU <- bU - b1
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
        cL <- abs((ml$maximum) - bL) > e
        cU <- abs((ml$maximum) - bU) > e
      }
      if (cU == "FALSE"){
        bL <- bL + b1
        bU <- bU + b1
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
        cL <- abs((ml$maximum) - bL) > e
        cU <- abs((ml$maximum) - bU) > e
      }
      while ("NaN" %in% ml$objective) {
        b1 <- b1*0.25
        bL <- -(b1) #lower bound (adjusted)
        bU <- b1  #upper bound (adjusted)
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
      }
    }
    
    B <- ml$maximum
    maxl <- ml$objective
    
    #Output results 
    out <- data.frame(chr, windowsize, window, s ,loc, maxl, B, hap_N, nM0, nF0, nM1, nF1, p, q, c)
    colnames(out) <- (c("Chromosome","Window_Size", "Window", "Site", "Locus", "ML", "Beta", "N_haplotypes", "N_males_0", "N_females_0", "N_males_1", "N_females_1","p","q","c"))
    output <- rbind(output, out)
    
  }
  return(output)
}

#optimization for bootstrapping (modify output)
likelihood.out.boot <- function(data,windowsize,window){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- unique(as.numeric(sub("[A-z].", "", data[data$Window == window,2])))
  
  #assign values
  for (s in 1:windowsize){    #site by site
    
    #variables
    frq <- NULL
    numM <- NULL
    numF <-NULL
    h <- NULL
    f <- NULL
    nM_0 <- NULL
    nF_0 <- NULL
    
    for (i in haplist){   
      frq[i+1] <- c(sex_avg_frq(data, window,i))
      numM[i+1] <- c(haplonum(data, i, window, "m"))
      numF[i+1] <- c(haplonum(data, i, window, "f"))
      
      if (frq[i+1] > 0){
        h[i+1] <- c(((numM[i+1]+numF[i+1])*log(frq[i+1])))
      }
      k = ((windowsize+1) - s)
      if (allele(i,k) == 0){  
        f[i+1] <- frq[i+1] 
        nM_0[i+1] <- numM[i+1] 
        nF_0[i+1] <- numF[i+1]
      }
      
    }
    
    #Calculate variables
    nM0 = sum(nM_0, na.rm=TRUE)
    nF0 = sum(nF_0, na.rm=TRUE)
    nM1 = (sum(numM, na.rm=TRUE) - sum(nM0, na.rm=TRUE))
    nF1 = (sum(numF, na.rm=TRUE) - sum(nF0, na.rm=TRUE))
    c = sum(h, na.rm=TRUE)
    loc = colnames(hap_test[wins[[window]][s]])
    q = ((nM0/(nM0+nM1)) + (nF0/(nF0+nF1)))/2
    p = 1-q
    
    #Optimize the likelihood
    hap_counts <- c(nF0,nM0,nF1,nM1,c)
    allele_freqs <-c(p,q)
    
    #Dynamic bounds for LL search
    b1 <- 0.1
    bL <- -(b1)   #lower bound
    bU <- b1    #upper bound
    e <- 10^(-4)
    
    ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
    
    while ("NaN" %in% ml$objective) {
      b1 <- b1*0.25
      bL <- -(b1)   #lower bound (adjusted)
      bU <- b1    #upper bound (adjusted)
      ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
    }
    
    cL <- abs((ml$maximum) - bL) > e
    cU <- abs((ml$maximum) - bU) > e
    
    while ( cL == "FALSE" | cU == "FALSE"){
      if (cL == "FALSE"){
        bL <- bL - b1
        bU <- bU - b1
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
        cL <- abs((ml$maximum) - bL) > e
        cU <- abs((ml$maximum) - bU) > e
      }
      if (cU == "FALSE"){
        bL <- bL + b1
        bU <- bU + b1
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
        cL <- abs((ml$maximum) - bL) > e
        cU <- abs((ml$maximum) - bU) > e
      }
      while ("NaN" %in% ml$objective) {
        b1 <- b1*0.25
        bL <- -(b1) #lower bound (adjusted)
        bU <- b1  #upper bound (adjusted)
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
      }
    }
    
    B <- ml$maximum
    maxl <- ml$objective
    
    #Output results 
    out <- data.frame(window,loc,maxl,B)
    colnames(out) <- (c("Window","Locus","ML","Boot_Beta"))
    output <- rbind(output, out)
    
  }
  return(output)
}

#Optimization for permuations
likelihood.out.rand <- function(data,windowsize,window){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- unique(as.numeric(sub("[A-z].", "", data[data$Window == window,2])))
  
  #assign values
  for (s in 1:windowsize){    #site by site
    
    #variables
    frq <- NULL
    numM <- NULL
    numF <-NULL
    h <- NULL
    f <- NULL
    nM_0 <- NULL
    nF_0 <- NULL
    
    for (i in haplist){   
      frq[i+1] <- c(sex_avg_frq(data, window,i))
      numM[i+1] <- c(haplonum(data, i, window, "m"))
      numF[i+1] <- c(haplonum(data, i, window, "f"))
      
      if (frq[i+1] > 0){
        h[i+1] <- c(((numM[i+1]+numF[i+1])*log(frq[i+1])))
      }
      k = ((windowsize+1) - s)
      if (allele(i,k) == 0){  
        f[i+1] <- frq[i+1] 
        nM_0[i+1] <- numM[i+1] 
        nF_0[i+1] <- numF[i+1]
      }
      
    }
    
    #Calculate variables
    nM0 = sum(nM_0, na.rm=TRUE)
    nF0 = sum(nF_0, na.rm=TRUE)
    nM1 = (sum(numM, na.rm=TRUE) - sum(nM0, na.rm=TRUE))
    nF1 = (sum(numF, na.rm=TRUE) - sum(nF0, na.rm=TRUE))
    c = sum(h, na.rm=TRUE)
    loc = colnames(hap_test[wins[[window]][s]])
    q = ((nM0/(nM0+nM1)) + (nF0/(nF0+nF1)))/2
    p = 1-q
    
    #Optimize the likelihood
    hap_counts <- c(nF0,nM0,nF1,nM1,c)
    allele_freqs <-c(p,q)
    
    #Dynamic bounds for LL search
    b1 <- 0.1
    bL <- -(b1)   #lower bound
    bU <- b1    #upper bound
    e <- 10^(-4)
    
    ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
    
    while ("NaN" %in% ml$objective) {
      b1 <- b1*0.25
      bL <- -(b1)   #lower bound (adjusted)
      bU <- b1    #upper bound (adjusted)
      ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
    }
    
    cL <- abs((ml$maximum) - bL) > e
    cU <- abs((ml$maximum) - bU) > e
    
    while ( cL == "FALSE" | cU == "FALSE"){
      if (cL == "FALSE"){
        bL <- bL - b1
        bU <- bU - b1
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
        cL <- abs((ml$maximum) - bL) > e
        cU <- abs((ml$maximum) - bU) > e
      }
      if (cU == "FALSE"){
        bL <- bL + b1
        bU <- bU + b1
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
        cL <- abs((ml$maximum) - bL) > e
        cU <- abs((ml$maximum) - bU) > e
      }
      while ("NaN" %in% ml$objective) {
        b1 <- b1*0.25
        bL <- -(b1) #lower bound (adjusted)
        bU <- b1  #upper bound (adjusted)
        ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
      }
    }
    
    B <- ml$maximum
    maxl <- ml$objective
    
    #Output results 
    out <- data.frame(window,loc,maxl,B)
    colnames(out) <- (c("Window","Locus","ML","Rand_Beta"))
    output <- rbind(output, out)
    
  }
  return(output)
}


#Take a bootstrap sample
bootstrap <- function(data){
  
  polycheck_m = FALSE   #check polymorphism
  polycheck_f = FALSE    
  
  while (polycheck_m == FALSE | polycheck_f == FALSE) {
    
    # Bootstrap haplotypes for males and females within each window
    males_boot <- lapply(males, sample, replace = TRUE) #1 for males
    females_boot <- lapply(females, sample, replace = TRUE) #2 for females
    
    polycheck_m <- all(polycheck_haps(males))
    polycheck_f <- all(polycheck_haps(females))
  }
  
  #Remove Ambiguous Sites (with ?)
  #males <- lapply(males, function(x) x[!grepl("?", x)])
  #females <- lapply(females, function(x) x[!grepl("?", x)])
  
  #Count the haplotypes and convert to base 10
  haplocounts(males_boot, 1) -> males_counts
  haplocounts(females_boot, 2) -> females_counts
  
  hapcounts <- data.frame(Window = 1:length(wins))
  colnamevec <- c("haplotype", "count")
  hapcounts[ ,colnamevec] <- NA
  
  for (i in 1:length(females_counts)) {
    females_counts[[i]]$Window <- as.numeric(names(females_counts[i]))
  }
  for (i in 1:length(males_counts)) {
    males_counts[[i]]$Window <- as.numeric(names(males_counts[i]))
  }
  
  female_haps <- do.call(rbind, females_counts)
  male_haps <- do.call(rbind, males_counts)
  female_haps$haplotype <- sub("^", "f.", female_haps$haplotype)
  male_haps$haplotype <- sub("^", "m.", male_haps$haplotype)
  
  #Create merged dataframe for haplotype counts
  merge(hapcounts, female_haps, by=c("Window","haplotype"), all=TRUE) -> merged_df
  merge(merged_df, male_haps, by=c("Window","haplotype"), all=TRUE) -> merged_df
  merged_df$count.y[is.na(merged_df$count.y)] <- merged_df$count[is.na(merged_df$count.y)]
  merged_df %>%
    select(-count.x,-count) -> haplotype_counts
  na.omit(haplotype_counts) -> haplotype_counts
  names(haplotype_counts)[names(haplotype_counts) == 'count.y'] <- 'counts'
  
  return(haplotype_counts)
}

#Perform likelihood calculations with a bootstrap sample and print to sdout
boot_data <- function(data){
  haplotype_counts_boot <- bootstrap(data)
  boot_betas <- NULL
  
  #Likelihood calculations
  winlist <- names(wins)
  boot_ml_output <- lapply(winlist, likelihood.out.boot, data = haplotype_counts_boot, windowsize = win)
  boot_ml_output <- do.call(rbind, boot_ml_output)
  boot_ml_output$Locus <- gsub('\\s+', '', boot_ml_output$Locus)
  boot_ml_output$Locus <- as.numeric(boot_ml_output$Locus)
  boot_ml_output[order(boot_ml_output$Window,boot_ml_output$Locus),] -> boot_ml_output
  boot_ml_output %>%
    select(Window,Locus,ML,Boot_Beta) -> boot_ml_output
  boot_betas <- cbind(boot_betas, boot_ml_output$Boot_Beta)
  boot_betas <- as.matrix(t(boot_betas))
  cat(boot_betas, "\n")
}

#permutation function for haplotype counts
permute_sample <- function(data){
  
  polycheck_m = FALSE   #check polymorphism
  polycheck_f = FALSE    
  
  while (polycheck_m == FALSE | polycheck_f == FALSE) {
    
    suppressMessages(perm_data <- data %>% 
                       distinct(sample, sex) %>% 
                       mutate(perm_sex = sample(sex)) %>% 
                       inner_join(hap_test) %>%
                       select(-sex, sex=perm_sex))
    
    data_m <- perm_data %>% filter(sex=="1") 
    data_f <- perm_data %>% filter(sex=="2")
    
    data_new <- rbind(data_m, data_f)
    
    polycheck_m <- all(apply(data_m[3:length(data_m)], 2, function(x) length(unique(x))) != 1)
    polycheck_f <- all(apply(data_f[3:length(data_f)], 2, function(x) length(unique(x))) != 1)
  }
  
  # Get haplotypes for males and females
  males <- lapply(wins,gethaps,input=data_new,s=1) #1 for males
  females <- lapply(wins,gethaps,input=data_new,s=2) #2 for females
  
  #Remove Ambiguous Sites (with ?)
  #males <- lapply(males, function(x) x[!grepl("?", x)])
  #females <- lapply(females, function(x) x[!grepl("?", x)])
  
  #Count the haplotypes and convert to base 10
  haplocounts(males, 1) -> males_counts
  haplocounts(females, 2) -> females_counts
  
  hapcounts <- data.frame(Window = 1:length(wins))
  colnamevec <- c("haplotype", "count")
  hapcounts[ ,colnamevec] <- NA
  
  for (i in 1:length(females_counts)) {
    females_counts[[i]]$Window <- as.numeric(names(females_counts[i]))
  }
  for (i in 1:length(males_counts)) {
    males_counts[[i]]$Window <- as.numeric(names(males_counts[i]))
  }
  
  female_haps <- do.call(rbind, females_counts)
  male_haps <- do.call(rbind, males_counts)
  female_haps$haplotype <- sub("^", "f.", female_haps$haplotype)
  male_haps$haplotype <- sub("^", "m.", male_haps$haplotype)
  
  #Create merged dataframe for haplotype counts
  merge(hapcounts, female_haps, by=c("Window","haplotype"), all=TRUE) -> merged_df
  merge(merged_df, male_haps, by=c("Window","haplotype"), all=TRUE) -> merged_df
  merged_df$count.y[is.na(merged_df$count.y)] <- merged_df$count[is.na(merged_df$count.y)]
  merged_df %>%
    select(-count.x,-count) -> haplotype_counts
  na.omit(haplotype_counts) -> haplotype_counts
  names(haplotype_counts)[names(haplotype_counts) == 'count.y'] <- 'counts'
  
  return(haplotype_counts)
}

#get beta from permuted sample
rand_data <- function(data){
  haplotype_counts_rand <- permute_sample(data)
  rand_betas <- NULL
  
  #Likelihood calculations
  winlist <- names(wins)
  rand_ml_output <- lapply(winlist, likelihood.out.rand, data = haplotype_counts_rand, windowsize = win)
  rand_ml_output <- do.call(rbind, rand_ml_output)
  rand_ml_output$Locus <- gsub('\\s+', '', rand_ml_output$Locus)
  rand_ml_output$Locus <- as.numeric(rand_ml_output$Locus)
  rand_ml_output[order(rand_ml_output$Window),] -> rand_ml_output
  rand_ml_output %>%
    select(Window,Locus,ML,Rand_Beta) -> rand_ml_output
  rand_betas <- cbind(rand_betas, rand_ml_output$Rand_Beta)
  rand_betas <- as.matrix(t(rand_betas))
  cat(rand_betas, "\n")
}