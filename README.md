
<!-- badges: start -->
<!-- badges: end -->

# poststratIV

Package to post-stratify an IV estimate to increase precision. Also has
replication code for associated paper, including code to run the
simulations and generate synthetic RCT with non-compliance data.

## Demo of use

First we load some data (or, in this case, make some with our built-in
DGP function):

``` r
library( poststratIV )
dat = make_dat( N = 1000,
                params = default_sim_params(),
                include_POs = FALSE )
head( dat )
#> # A tibble: 6 × 5
#>   X     complier     Z     S  Yobs
#>   <chr>    <dbl> <dbl> <dbl> <dbl>
#> 1 X2           0     0     0  2.11
#> 2 X2           0     0     0  3.13
#> 3 X2           0     0     0  3.40
#> 4 X4           1     0     0  2.26
#> 5 X1           0     1     0  3.88
#> 6 X1           0     1     0  4.80
```

The data needs a single categorical covariate to stratify on.

Once we have our data, we can estimate, declaring the outcome,
compliance indicator, treatment assignment indicator, and stratification
variable (baseline categorical variable of some sort):

``` r
res <- IV.est.strat( dat, 
                     Yobs="Yobs", S="S", Z="Z", strat_var="X" )
res %>%
    dplyr::select( Xblk, ITT.hat, pi.hat, LATE.hat, SE.wald ) %>%
    knitr::kable( digits = 2 )
```

| Xblk    | ITT.hat | pi.hat | LATE.hat | SE.wald |
|:--------|--------:|-------:|---------:|--------:|
| X1      |   -0.10 |   0.00 |     0.00 |    0.00 |
| X2      |    0.05 |   0.05 |     0.98 |    1.37 |
| X3      |    0.62 |   0.64 |     0.96 |    0.25 |
| X4      |    1.50 |   0.95 |     1.58 |    0.24 |
| UNSTRAT |    0.07 |   0.16 |     0.45 |    0.39 |
| IV_w    |    0.28 |   0.23 |     1.20 |    0.26 |
| DSS0    |    0.28 |   0.23 |     1.20 |    0.26 |
| IV_a    |    0.12 |   0.14 |     0.90 |    0.32 |
| DSS     |    0.28 |   0.23 |     1.20 |    0.26 |
| PWIV    |    0.28 |   0.23 |     1.28 |    0.17 |
| DSF     |    0.28 |   0.23 |     1.20 |    0.26 |

There is also the usual IV with no covariate adjustment:

``` r
res <- IV.est( dat, 
                     Yobs="Yobs", S="S", Z="Z" )
knitr::kable( res, digits=2 )
```

| ITT.hat | pi.hat | LATE.hat | SE.ITT | SE.pi | cov.ITT.pi | SE.wald | SE.delta |    n |
|--------:|-------:|---------:|-------:|------:|-----------:|--------:|---------:|-----:|
|    0.07 |   0.16 |     0.45 |   0.06 |  0.02 |          0 |    0.39 |     0.41 | 1000 |

## Replication Simulation Code

The `replication` directory holds all the replication files for the
post-stratification IV paper.

To replicate the results, first install the poststratIV package, which
has all the code for the estimators and the DGP.

    devtools::install_github( "lmiratrix/poststratIV" )

You then download the `replication` folder from GitHub. You can download
just the replication folder via [the
gitzip](https://kinolien.github.io/gitzip/) web interface, by navigating
to it and then putting the following in the search bar:

    https://github.com/lmiratrix/poststratIV/tree/main/replication

Once you download the replication folder, open the R project inside the
folder (so the `here::here` commands work correctly) and then you need
to run the scripts in the folder in the following order:

- The simulation files (any order). Note that you will run both the main
  simulation and the predC simulation scripts twice, once with
  ONE_SIDED_SIMULATION=TRUE and once with it to FALSE. This runs
  everything with one- or two-sided noncompliance. It saves the results
  in two different directories. The analysis files will then stack these
  results to give a unified analysis.

- The analysis files (any order)

The scripts will build directories as needed, saving results in those
directories.

You can run all of it via the `run_everything.R` script.

## Replication Applied Example Code

This code is also in `replication`. You will need to download the data
and put it in a `data` directory in the `replication` folder. Data can
be found at (<https://isps.yale.edu/research/data/d017>).
