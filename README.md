LQAS-IMP
==========

This repository provides code that accompanies the *Adapting Lot Quality Assurance Sampling (LQAS) to accommodate imperfect tests: application to COVID-19 serosurveillance in Haiti* paper. 

Setup
------------
The goal is to classify an area as “high” or “low” on some trait of interest. In our paper, we were interested in classifying hospitals as having high or low COVID-19 antibody prevalence among healthcare workers in Haiti. In the traditional form of LQAS, the sample size, _n_, and decision rule, _d_, are determined by the population size, _N_, and four parameters defined by users based on the specific context and goals. 

The four parameters are: 
+ *p.upper* the upper prevalence threshold by which an area is classified as high
+ *p.lower* the lower prevalence threshold by which an area is classified as low
+ *alpha.constraint* the probability that a high area is mistakenly classified as low
+ *beta.constraint* the probability that a low area is mistakenly classified as high

As our procedure also accounts for imperfect testing properties, we require that the *sensitivity* and *specificity* of the test is known to appropriately determine the LQAS system (_n_ and _d_).

Example
------------

Download the GitHub repository and source the R code.  
``` r
source("R/functions.R")
```

Use the `lqas_imp` function to determine the required sample size and decision rule:  
``` r
lqas_imp(N=80,
         p.upper=0.40,
         p.lower=0.10,
         sens=0.98,
         spec=0.95,
         alpha.contraint=0.15,
         beta.constraint=0.05
         )
```

The required sample size is 17 and the decision rule is 5. If at least 5 health workers test positive for COVID-19, then we would classify the facility as high prevalence. 

Computing time
------------
The above function will take considerably longer to run as the population size N increases. For example, the above code took about 3 minutes and 33 seconds to run. For any population size above 100, we recommend running the function in parallel. Specifically, one can calculate the required sample size (n) without adjustment for sensitivity and specifcity using the `lqas_hgm` function here or the many online tools available. This will provide a lower bound for the search space for the sample size n. Then, the function can be run in parallel for small increments of the sample size n and decision rule d. 

We are currently creating an R package that will reduce computing time by conducting a grid search. Please stay tuned. 
