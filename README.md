
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
#> 1 X1           0     0     0 4.16 
#> 2 X3           1     0     0 1.99 
#> 3 X1           0     1     0 5.23 
#> 4 X3           0     1     0 2.15 
#> 5 X1           0     0     0 3.40 
#> 6 X4           1     0     0 0.841
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
| X1      |    0.01 |   0.01 |     1.78 |    8.53 |
| X2      |   -0.10 |   0.05 |    -2.09 |    1.46 |
| X3      |    0.77 |   0.52 |     1.46 |    0.39 |
| X4      |    1.73 |   1.00 |     1.73 |    0.15 |
| UNSTRAT |    0.18 |   0.11 |     1.56 |    0.56 |
| IV_w    |    0.12 |   0.12 |     0.99 |    0.36 |
| DSS0    |    0.12 |   0.12 |     0.99 |    0.36 |
| IV_a    |    0.12 |   0.12 |     0.99 |    0.36 |
| DSS     |    0.20 |   0.21 |     0.97 |    0.30 |
| PWIV    |    0.12 |   0.12 |     1.66 |    0.14 |
| DSF     |    0.20 |   0.21 |     0.97 |    0.30 |

There is also the usual IV with no covariate adjustment:

``` r
res <- IV.est( dat, 
                     Yobs="Yobs", S="S", Z="Z" )
knitr::kable( res, digits=2 )
```

| ITT.hat | pi.hat | LATE.hat | SE.ITT | SE.pi | cov.ITT.pi | SE.wald | SE.delta |    n |
|--------:|-------:|---------:|-------:|------:|-----------:|--------:|---------:|-----:|
|    0.18 |   0.11 |     1.56 |   0.06 |  0.02 |          0 |    0.56 |     0.63 | 1000 |

## Replication Code

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
