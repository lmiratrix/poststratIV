

## This script focuses a simulation on a specific context
##
## Useful for more targeted explorations, etc.
##
## Script adapted from simulation_study_v2.R
##
## It can look at one- or two-sided noncompliance by setting the flag
## below.

library( tidyverse )
library( poststratIV )

# Run in parallel?
PARALLEL = TRUE


# Set TRUE for one-sided simulation
if ( !exists( "ONE_SIDED_SIMULATION" ) ) {
    ONE_SIDED_SIMULATION = TRUE
}

library( cli )

cli_alert_info( "Running the predictive-of-compliance simulation with ONE_SIDED_SIMULATION={ONE_SIDED_SIMULATION}\n" )



if ( !exists( "R" ) ) {
    # Number of sim runs
    R = 6 # 3000

    # chunkNo is a hack to make a bunch of smaller chunks for doing parallel more
    # efficiently.
    M_CHUNK = 2 #5
}


if ( ONE_SIDED_SIMULATION ) {
    result_dir = here::here( "results_onesided/" )
} else {
    result_dir = here::here( "results_twosided/" )
}

if ( !dir.exists( result_dir ) ) {
    dir.create( result_dir )
}


# Make directory for where to store results
make_result_dir_name <- function( name ) {
    paste0( result_dir, name )
}





factors = expand_grid( chunkNo = 1:M_CHUNK,
                       N = c( 3000, 10000 ),
                       pi_c = c( 0.15 ),
                       nt_shift = c( 0 ),
                       pred_comp = c(  "yes" ),
                       pred_Y = c( "no" ),
                       het_tx = c(  "no", "yes" ),
                       one_sided = ONE_SIDED_SIMULATION,
                       perfect_X = FALSE,
                       scaled_C = seq( 0, 1, length.out=7 )^3,
                       sd0 = 1 )
factors
factors <- factors %>% mutate(
    reps = R / M_CHUNK,
    seed = 36200329 + 1:n()
)
factors$chunkNo = NULL

nrow( factors )

pstrat = c( 0.35, 0.35, 0.25, 0.15 )
tau = c( 0, 0.2, 0.4, 0.6 )


# Print out what we are dealing with.
ff = factors[27,]
ff
ll = as.list( ff )
ll$reps = ll$seed = ll$chunkNo = NULL
ll$pstrat = pstrat
ll
#ll$N = 100000
dd = do.call( make_dat, ll )
cat( "\nData description of example dataset:\n" )
print( describe_sim_data(dd) )
summarize_sim_data(dd) %>%
    print()

dd$grp = paste0( dd$Z, "-", dd$complier )
ggplot( dd, aes( X, Yobs, col=grp ) ) +
    geom_boxplot()
sd( dd$Yobs )
mean( dd$Yobs )

if ( FALSE ) {
    ss <- run_sim( N = 1000, pi_c = 0.15, scaled_C = 1, pstrat = pstrat, tau = tau, sd0 = 1,
                   nt_shift = 0, params = NULL,
                   reps=10, pred_comp="yes", pred_Y ="no", het_tx="yes" )
    ss
    ss %>%
        group_by( method ) %>%
        summarise( bias = mean( LATE.hat - LATE, na.rm=TRUE ) )
}



#### Run the simulation #####


# Do the run!

tictoc::tic()

# parallel version
if ( PARALLEL ) {


    parallel::detectCores()

    library(future)
    library(furrr)

    n_worker = parallel::detectCores() - 2
    cat( "Using", n_worker, "workers\n" )

    plan(multisession, workers = n_worker )
    safe_run_sim = safely( run_sim )

    factors$res <- future_pmap( factors, .f = safe_run_sim,
                                pstrat = pstrat, tau = tau, params = NULL,
                                .options = furrr_options(seed = NULL),
                                .progress = TRUE )
} else {
    factors$res <- pmap( factors, .f = safe_run_sim,
                         pstrat = pstrat, tau = tau, params = NULL,
                         .progress = TRUE )

}
tictoc::toc()

factors
head( factors$res[[1]] )

sim_results <-
    factors %>%
    unnest(cols = res)


sim_results$sr = rep( c("res","err"), nrow(sim_results)/2)
sim_results = pivot_wider( sim_results, names_from = sr, values_from = res )

a = sim_results$err
lgc = map_lgl( a, is.null )
table( lgc )
lgc = map_dbl( a, length )
table( lgc )

print( sim_results )

cat( "First error message as sanity check:\n" )
print( sim_results$err[[1]] )

cat( "Saving results\n" )

saveRDS( sim_results,
         file=make_result_dir_name( "simulation_results_predC.rds" ) )

cat( "pred_C Script finished\n" )


