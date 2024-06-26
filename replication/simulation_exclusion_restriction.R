

## This script looks at how violations of the exclusion restriction is
## mitigated if we can pack compilers into a single strata.

library( tidyverse )
library( poststratIV )
library( cli )

# Run in parallel?
PARALLEL = TRUE

cat( "Running the exclusion restriction simulation\n" )


if ( !exists( "R" ) ) {
    # Number of sim runs
    warning( "Setting R=6 to demo code.  Raise number of simulation replicates for real results" )
    R = 6

}

result_dir = here::here( "results/" )
if ( !dir.exists( result_dir ) ) {
    dir.create( result_dir )
}

# Make directory for where to store results
make_result_dir_name <- function( name ) {
    paste0( result_dir, name )
}


cat( "making simulation scenarios\n" )


# chunkNo is a hack to make a bunch of smaller chunks for doing parallel more
# efficiently.
M_CHUNK = 10
factors = expand_grid( chunkNo = 1:M_CHUNK,
                       N = c( 3000 ),
                       pi_c = c( 0.15 ),
                       nt_shift = c( 0 ),
                       pred_comp = c(  "yes" ),
                       pred_Y = c( "no" ),
                       het_tx = c(  "no" ),
                       one_sided = FALSE,
                       exclusion_impact = 0.2,
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

# Set up strata sizes and tau
pstrat = c( 0.35, 0.35, 0.25, 0.15 )
tau = c( 0.5, 0.5, 0.5, 0.5 ) #, 0.2, 0.4, 0.6 )



# Print out what we are dealing with.
ff = factors[1,]
ff
ll = as.list( ff )
ll$reps = ll$seed = ll$chunkNo = NULL
ll$pstrat = pstrat
ll$tau = tau
ll
#ll$N = 100000
dd = do.call( make_dat, ll )
cat( "\nData description of example dataset:\n" )
print( describe_sim_data(dd) )
poststratIV:::summarize_sim_data(dd) %>%
    print()



#### Testing ####

if ( FALSE ) {

    f = factors[5,]

    rrr <- pmap( f, .f = run_sim,
                         pstrat = pstrat, tau = tau, params = NULL,
                         .progress = TRUE )
    rrr = rrr[[1]]
    summary( rrr$LATE )
    summary( rrr$ITT )
    0.5 * 0.15 + 0.2 * (1-0.15)

}


#### Run the simulation #####


# Do the run!

tictoc::tic()


cat( "There are", nrow(factors), "chunks of", factors$reps[[1]], "iterations each\n" )

safe_run_sim = safely( run_sim )

if ( PARALLEL ) {
    # parallel version

    parallel::detectCores()

    library(future)
    library(furrr)

    n_worker = parallel::detectCores() - 2
    cat( "Using", n_worker, "workers\n" )

    plan(multisession, workers = n_worker )

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
         file=make_result_dir_name( "simulation_results_exclusion.rds" ) )



cat( "exclusion restriction Script finished\n" )


