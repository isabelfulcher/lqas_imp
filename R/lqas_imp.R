library(data.table)
library(dplyr)
library(pbapply)


lqas_imp <- function(N, # population size
                     p.upper, # designate as "high" prevalence
                     p.lower, # designate as "low" prevalence
                     sens, # sensitivity of the test
                     spec, # specifciity of the test
                     alpha.contraint, # alpha 
                     beta.constraint, # beta
                     ){
  
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


