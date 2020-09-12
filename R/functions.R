library(data.table)
library(dplyr)
library(pbapply)

# LQAS classification system for imperfect tests
lqas_imp <- function(N, 
                     p.upper, 
                     p.lower, 
                     sens, 
                     spec, 
                     alpha.contraint,
                     beta.constraint){
  
  N.up.pos<-round(p.upper*N)
  N.up.neg<-N-N.up.pos
  N.low.pos<-round(p.lower*N)
  N.low.neg<-N-N.low.pos
  
  compute_vals <- function(d, n){
    
    dt <- data.table(expand.grid(i=0:n, k=0:n, l=0:n))
    dt <- dt[l<=k]
    dt[, tp := dbinom(l,size=k,sens)]
    dt[, fp := dbinom(i-l,size=n-k,1-spec)]
    dt[, pl := dhyper(x=k,m=N.low.pos,n=N.low.neg,k=n)]
    dt[, pu := dhyper(x=k,m=N.up.pos,n=N.up.neg,k=n)]
    dt[,product_low := tp*fp*pl]
    dt[,product_up := tp*fp*pu]
    
    data.frame(
      n,
      d,
      alpha = sum(dt[i < d][["product_up"]]),
      beta = sum(dt[i >= d][["product_low"]])
    )
  }
  
  
  # Populate possible values and then evaluate them systematically
  dn <- expand.grid(d=1:N,n=1:N) %>% dplyr::filter(d <= n)
  set.seed(1)
  pblapply(sample(1:nrow(dn)), function(idx){
    compute_vals(dn[idx,"d"], dn[idx,"n"])
  }) %>% rbindlist() %>% data.frame() -> value_df
  
  choice <- value_df %>% filter(alpha <= 0.10 & beta <=.10) %>% filter(n==min(n)) 
  n_min <- choice$n
  d_min <- choice$d
  
  # return LQAS system
  return(data.frame(N=N,n=n_min,d=d_min))
  
}


# Classic LQAS classification system (Code written by Bethany Hedt-Gauthier)
lqas_hgm <- function(p_upper,p_lower,alpha.constraint,beta.constraint,N,sens,spec){
  
  #adjust for sens and spec of test
  p.low_adj<-p_lower*sens+(1-p_lower)*(1-spec)
  p.up_adj<-p_upper*sens+(1-p_upper)*(1-spec)
  
  #index possible number of positive tests 
  n<-c(1:N)
  
  #number of successes and failures in finite population of size N at upper and lower threshold 
  ## need for hypergeomtetric distribution input
  N.up.pos<-round(p.up_adj*N)
  N.up.neg<-N-N.up.pos
  N.low.pos<-round(p.low_adj*N)
  N.low.neg<-N-N.low.pos
  
  # Storage vector, populate with -5
  ## first row is potential sample size, n
  ## second row is potential decision rule (indexed by x below)
  results<-matrix(-5,nrow=2, ncol=length(n)) 
  results[1,]<-n
  
  # loop over all possible sample sizes and find potential decision rules
  for(j in 1:length(n)){
    
    #possible number of cases
    x<-c(0:n[j])
    
    #initialize alpha and beta vectors
    alpha<-rep(-5,length(x)+1)
    beta<-rep(-5,length(x)+1)
    
    alpha[1]<-0
    beta[1]<-1
    
    alpha[length(x)+1]<-1
    beta[length(x)+1]<-0
    
    #For sample size k, evaluate probability of observing x for p_lower and p_upper
    ## use dhyper (x,m,n,k) where x is number of successes because FINITE population of size N
    ## m is total number of successes in population
    ## n is total number of failures in population
    ## k is sample size (number pulling from population)
    alpha.prel<-dhyper(x,m=N.up.pos,n=N.up.neg,k=n[j]) #if you choose k under the UPPER limit, how likely are you to see x cases?
    beta.prel<-dhyper(x,m=N.low.pos,n=N.low.neg,k=n[j]) #if you choose k under the LOWER limit, how likely are you to see x cases?
    
    # store appropriate alpha and beta values at each x with sample size k
    for(i in 2:length(x)){
      alpha[i]<-sum(alpha.prel[1:(i-1)]) #P(X < d | p = pu)
      beta[i]<-sum(beta.prel[i:length(x)]) #P(X >= d | p = pl)
    }
    
    # which values of x meet alpha and beta constraints
    condition1<-which(alpha<alpha.constraint)
    condition2<-which(beta<beta.constraint)
    
    # if both conditions are met then enter value of x (decision rule) --> candidate!
    if(max(condition1)==min(condition2)){
      results[2,j]<-max(condition1)-1
    }
  }
  
  # Return all candidate sample sizes with decision rules
  possibles<-which(results[2,]>-1)
  
  # Choose minimum sample size and corresponding decision rule
  min.n<-results[1,min(possibles)]
  cor.d<-results[2,min(possibles)]
  
  # Return finite pop sample, minimum sample size, and decision rule
  return(c(N,min.n,cor.d))
}


