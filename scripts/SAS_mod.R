#simSAS-modified
simSAS <- function(p, f, n, k, beta, samp_m, samp_f, num) {
  #went and corrected samp to be an equal 50-50 split 
  #p = minor allele freq
  #r^2 = corr between sites
  #f = frequency of each haplotype of length n
  #k = site under selection, site 3
  #n = number of sites in haplotype
  #beta = strength of selection
  #samp_m = men in sample
  #samp_f = females in sample
  #num = number of samples, this can be set to 1 
  #simSAS_start = Sys.time()
  haplos <- getHaplos(n)
  # Make list of frequencies for all haplotypes possible of length n
  freqs <- f
  #Getting male and female haplotype frequencies after selection for all possible haplotypes
  fhm  <- sapply(haplos, function(x) getMaleHap(freqs[x + 1], x, beta, k, p))
  fhf  <- sapply(haplos, function(x) getFemaleHap(freqs[x + 1], x, beta, k, p))
  
  #Getting data that reflects those haplotype frequencies sampled from multinomial distribution.
  lnmales <- rmultinom(num, size = samp_m, prob = fhm)
  #lnmales <- rmultinom(num, size = (samp/2), prob = fhm)
  
  lnfemales <- rmultinom(num, size = samp_f, prob = fhf)
  #lnfemales <- rmultinom(num, size = (samp/2), prob = fhf)
  
  haps <- data_frame(haplos)
  fems <- data_frame(freqsF = fhf, lnfemales)
  mals <- data_frame(freqsM = fhm, lnmales)
  
  res <- c(haps, fems, mals)
  #SimSAS_end = Sys.time()
  #print(paste("SimSAS, 1 window:", SimSAS_end - simSAS_start))
  return(res)
}

getHaplos <- function(n){
  haplos <- c(0:((2^n) - 1))
  return(haplos)
}
getMaleHap <- function(fh, hap, beta, k, p){
  q <- 1-p
  if (getAllele(hap, k) == 0){
    fhm <- (1 + (p * beta)) * fh
  }else{
    fhm <- (1 - (q * beta)) * fh
  }
  return(fhm)
}
getAllele <- function(haplo, pos){
  allele <- floor((haplo/(2^(pos-1))) %% 2)
  return(allele)
}

getFemaleHap <- function(fh, hap, beta, k, p){
  q <- 1-p
  if (getAllele(hap, k) == 0){
    fhm <- (1 - (p * beta)) * fh
  }else{
    fhm <- (1 + (q * beta)) * fh
  }
  return(fhm)
}
