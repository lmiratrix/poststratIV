
<!-- badges: start -->
<!-- badges: end -->

# poststratIV

Package to post-stratify an IV estimate to increase precision. Also has
replication code for associated paper, including code to run the
simulations and generate synthetic RCT with non-compliance data.

## Demo of use

First we make or load some data:

``` r
library( poststratIV )
dat = make_dat( N = 1000,
                params = default_sim_params() )
```

Then we estimate, declaring the outcome, compliance indicator, treatment
assignment indicator, and stratification variable (baseline categorical
variable of some sort):

``` r
res <- IV.est.strat( dat, 
                     Yobs="Yobs", S="S", Z="Z", strat_var="X" )
res %>%
    dplyr::select( Xblk, ITT.hat, pi.hat, LATE.hat, SE.wald ) %>%
    knitr::kable( digits = 2 )
```

| Xblk    | ITT.hat | pi.hat | LATE.hat | SE.wald |
|:--------|--------:|-------:|---------:|--------:|
| X1      |    0.02 |   0.03 |     0.65 |    2.29 |
| X2      |    0.05 |   0.07 |     0.69 |    0.92 |
| X3      |    0.27 |   0.43 |     0.62 |    0.35 |
| X4      |    1.37 |   1.00 |     1.37 |    0.16 |
| UNSTRAT |    0.08 |   0.15 |     0.54 |    0.42 |
| IV_w    |    0.13 |   0.14 |     0.94 |    0.30 |
| DSS0    |    0.13 |   0.14 |     0.94 |    0.30 |
| IV_a    |    0.13 |   0.14 |     0.94 |    0.30 |
| DSS     |    0.13 |   0.14 |     0.94 |    0.30 |
| PWIV    |    0.13 |   0.14 |     1.22 |    0.15 |
| DSF     |    0.21 |   0.21 |     0.96 |    0.25 |
