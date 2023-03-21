#Optimization function site-by-site, vectorized 
#Updated: 12/7/2022 - JMC

#function
likelihood.out.vec <- function(data,windowsize,window){
  
  #total_lik_start = Sys.time()
  
  #output DF
  output <- data.table()
  
  #list of haplos from data
  haplist <- as.numeric(unique(unlist(data[data$Window == window,3])))
  
  #variables (pre-allocated)
  frq <-  vector("numeric", length(haplist))
  numM <- vector("numeric", length(haplist))
  numF <- vector("numeric", length(haplist))
  h <- vector("numeric", length(haplist))
  f <- vector("numeric", length(haplist))
  nM_0 <- vector("numeric", length(haplist))
  nF_0 <- vector("numeric", length(haplist))
  
  #get counts
  numM <- sapply(haplist, function(i) ifelse(length(data[data$Window == window & data$haplotype == i, ]$male_counts) == 0,0, data[data$Window == window & data$haplotype == i, ]$male_counts))
  numF <- sapply(haplist, function(i) ifelse(length(data[data$Window == window & data$haplotype == i, ]$female_counts) == 0,0, data[data$Window == window & data$haplotype == i, ]$female_counts))
  
  
  #allele counts
  k <- (windowsize + 1) / 2
  nM_0 <- sapply(haplist, function(i) ifelse(allele(i,k)==0, numM[i+1], 0))
  nF_0 <- sapply(haplist, function(i) ifelse(allele(i,k)==0, numF[i+1], 0))     
  nM0 = sum(nM_0, na.rm=TRUE)
  nF0 = sum(nF_0, na.rm=TRUE)
  nM1 = (sum(numM, na.rm=TRUE) - sum(nM0, na.rm=TRUE))
  nF1 = (sum(numF, na.rm=TRUE) - sum(nF0, na.rm=TRUE))
  
  
  #Calculate input variables and allele freqs
  #sex-avged freqs 
  frq <- sapply(haplist, function(i) {
    nMh = data[data$Window == window & data$haplotype == i, ]$male_counts
    nFh = data[data$Window == window & data$haplotype == i, ]$female_counts
    fM = nMh / (nM0+nM1)
    fF = nFh / (nF0+nF1)
    fq = (fM+fF)/2
    return(fq)
  })
  
  #null
  c = sum(((numM + numF) * log(frq)), na.rm=TRUE)
  
  qm = (nM0/(nM0+nM1))
  qf = (nF0/(nF0+nF1))
  q = (qm+qf)/2
  p = 1-q
  
  #Optimize the likelihood
  hap_counts <- c(nF0,nM0,nF1,nM1,c)
  allele_freqs <-c(p,q)
  
  #Dynamic bounds for LL search
  b1 <- 0.5
  bL <- -(b1)   #lower bound
  bU <- b1    #upper bound
  e <- 10^(-4)
  
  #opt_start_time = Sys.time()
  
  ml <- optimize(lik.function, lower=bL, upper=bU, maximum=TRUE, params = hap_counts, freqs = allele_freqs)
  
  #opt_end_time = Sys.time()
  
  #print(paste("Optimization Time for window", window, "is:", opt_end_time - opt_start_time))
  
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
  
  #get betas
  B<- ml$maximum
  maxl<- ml$objective
  
  #calculate p-values (LRT)
  LRT_stat<- -2 * (as.numeric(maxl) - as.numeric(c))
  p_val <- pchisq(abs(LRT_stat), df=1, lower.tail = FALSE)
  
  #output info
  chr <- unique(data$Chrom)
  pos <- unique(data[data$Window == window, ]$Center_Pos)
  snpid <- unique(data[data$Window == window, ]$Center_SNP)
  maf <- unique(data[data$Window == window, ]$Center_MAF)
  al0 <- unique(data[data$Window == window, ]$Center_Allele0)
  al1 <- unique(data[data$Window == window, ]$Center_Allele1)
  
  #Output results 
  out <- data.table(chr,window,k,pos,snpid,maf,al0,al1,nF0,nM0,nF1,nM1,p,q,c,maxl,B,p_val)
  colnames(out) <- (c("Chrom","Window","Site","Position","SNP","MAF","Allele0","Allele1","nF0","nM0","nF1","nM1","p","q","nullLL","maxLL","Beta_hat","p_value"))
  output <- rbind(output, out)
  #total_lik_end = Sys.time()
  #print(paste("Total Likelihood Time for window", window, "is:", total_lik_end - total_lik_start))
  return(output)
}