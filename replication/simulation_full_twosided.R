
# This simulation study, the main simulation of the paper, focuses on
# the LATE given a DGP that has CACE vary by proportion of complier.
#
# It examines two things:
#
# (1) The performance of the point estimate.
#
# (2) The performance of the standard error estimate
#     - do we get good coverage?
#     - Is the estimate of the SE estimate more or less unstable?
#
# This script saves the simulations as it goes along in a "frags"
# directory.  If the process dies, and then you re-run it, it will
# skip those parts of the simulation that are found in the directory.
# If you want to fully re-run the simulation, you will need to delete
# the "frags" directory entirely.



#### Script settings ####


# Set TRUE for one-sided simulation, FALSE for two-sided. Results will
# be saved in different directories, depending.
if ( !exists( "ONE_SIDED_SIMULATION" ) ) {
    ONE_SIDED_SIMULATION = TRUE
}


if ( !exists( "R" ) ) {
    # Number of sim runs
    R = 6
    warning( "Setting R=6 to demo code.  Raise number of simulation replicates for real results" )
    
    # M_CHUNK is a hack to make a bunch of smaller chunks for doing parallel more
    # efficiently.
    M_CHUNK = 2
}


# Run simulation in parallel?
PARALLEL = TRUE



#### Load libraries and set up directory names ####


library( tidyverse )
library( poststratIV )
library( cli ) # for printing to console


cli_alert_info( "Running the main simulation with ONE_SIDED_SIMULATION={ONE_SIDED_SIMULATION}\n" )


result_dir = NA
if ( ONE_SIDED_SIMULATION ) {
    result_dir = here::here( "results_onesided/" )
} else {
    result_dir = here::here( "results_twosided/" )
}


# Make directory for where to store results
make_result_dir_name <- function( name ) {
    paste0( result_dir, name )
}

if ( !file.exists( result_dir ) ) {
    dir.create(result_dir)
}

# Set up place to store partially run simulation results
frag_dir =  make_result_dir_name("frags/" )
if ( !file.exists( frag_dir ) ) {
    dir.create(frag_dir)
}






#### Run the multifactor simulation ####

# M_CHUNK is a hack to make a bunch of duplicate tasks for doing
# parallel more efficiently.  (I.e., we run the same scenario multiple
# times.)
factors = expand_grid( chunkNo = 1:M_CHUNK,
                       N = c( 500, 1000, 2000 ),
                       pi_c = c( 0.05, 0.075, 0.10 ),
                       nt_shift = c( -1, 0, 1 ),
                       pred_comp = c( "yes", "no" ),
                       pred_Y = c( "yes", "no" ),
                       het_tx = c( "yes", "no" ),
                       sd0 = 1, # c( 0.5, 1, 2 ),
                       one_sided = ONE_SIDED_SIMULATION
                       )
factors <- factors %>% mutate(
    reps = R / M_CHUNK,
    seed = 16200329 + 1:n()
)

nrow( factors )


if ( FALSE ) {

    f2 <- filter( factors,  N == 500, het_tx == "no", nt_shift == 0,
            pi_c == 0.05, pred_comp == "yes",  pred_Y=="no")
    f2
}


##### Debugging: Explore the scenarios used in the simulation #####
if ( FALSE ) {

    # Get scale of outcomes by generating some large sample datasets and
    # describing them, for all factor combinations.

    f2 = filter( factors, N == 2000 )

    f2 <- f2 %>%
        dplyr::select( -chunkNo, -seed, -reps ) %>%
        unique()
    nrow( f2 )

    # Make dataset of all scenarios
    f2$dat = pmap( f2, one_run, data_only=TRUE )
    f2 = rename( f2, sd0_param = sd0 )


    describe_sim_data( f2$dat[[1]] )
    summarize_sim_data( f2$dat[[1]]  )

    # Summarize the data
    f2$sum = map( f2$dat, describe_sim_data )
    f2$sum2 = map( f2$dat, summarize_sim_data )
    f2 = dplyr::select( f2, -dat )
    f2 = unnest( f2, cols=sum2 )
    f2 = unnest( f2, cols= sum,
                 names_repair = "unique")
    filter( f2, pred_comp == "yes", X != "ALL" ) %>%
        group_by( X ) %>%
        summarise( pi_mn = mean( pi ),
                   pi_q5 = quantile( pi, 0.05 ),
                   pi_q95 = quantile( pi, 0.95 ) )
    pull( pi ) %>% skimr::skim()
    ftop
    f2 = filter( f2, X == "ALL" ) %>%
        dplyr::select( -X, -n, -sd0_param, -N, -nt, -ITT, -per, -per_c  )
    f2

    # Look at explanatory power of covariates
    qplot( f2$R2_comp )
    summary( lm( R2_comp ~ pred_comp * (pred_Y + het_tx + pi_c), data=f2  ) )
    f2 %>% group_by( pred_comp, pi_c ) %>%
        summarise( sd_comp = sd( R2_comp ),
                   R2_comp = mean( R2_comp))

    qplot( f2$R2_y )
    summary( lm( R2_y ~ pred_Y * (pred_comp + het_tx + pi_c + nt_shift), data=f2  ) )

    f2 %>% group_by( pred_Y ) %>%
        summarise( R2_y = mean( R2_y))
    # Look at a few scenarios to see how things vary
    sample_n( f2, 10 )

    ### Looking at LATEs ###

    # Big LATE for het tx and complier prediction, due to necessity of
    # that.
    skimr::skim( f2$LATE )
    sample_n( f2, 16 ) %>%
        arrange( LATE )
    f2 %>% group_by( pred_comp, het_tx ) %>%
        summarise( mnLATE = mean( LATE ),
                   sdLATE = round( sd( LATE ), digits=4 )  )

    #filter( f2, pi_c == 0.10 )

    ### Looking at Y0 variance ###
    # Most of the control-side variances are about 0
    skimr::skim(f2$sd0)
    round( sort( f2$sd0 ), digits=2 )

    summary( lm( sd0 ~ (nt_shift + pi_c + pred_comp + pred_Y)^2, data=f2 ) )
}


###### Debugging: Run a single simulation trial #####

if ( FALSE ) {

    # Run a single simulation scenario, for debugging purposes
    ff = factors[5,]
    ff
    ll = as.list( ff )
    #ll$N = 100000
    dd = do.call( run_sim, ll )

    dd

    dd %>% group_by( method ) %>%
        summarise( meanNA = mean( is.na( LATE.hat ) ),
                   SE = sd( LATE.hat, na.rm=TRUE ) )

}





#### Do the full simulation run!   ####

cli_alert_info( "Running simulation in {nrow(factors)} chunks\n" )

tictoc::tic()


safe_run_sim = safely( run_sim )
file_saving_sim = function( chunkNo, seed, ... ) {
    fname = paste0( frag_dir, "fragment_", chunkNo, "_", seed, ".rds" )
    res = NA
    if ( !file.exists(fname) ) {
        res = safe_run_sim( seed=seed, ... )
        saveRDS(res, file = fname )
    } else {
        #scat( "\nSkipping chunk %s-%s\n", seed, chunkNo )
        res = readRDS( file=fname )
    }
    return( res )
}

# Shuffle the rows so we run the tasks in a random order to load
# balance (in case some scenarios are harder than others).
factors = sample_n( factors, nrow(factors) )

if ( FALSE ) {
    factors = filter( factors, N==500, pi_c == 0.05 )
    warning( "Dropped all non-small simulations" )
}

if ( PARALLEL ) {
    # Run in parallel

    parallel::detectCores()

    library(future)
    library(furrr)

    #plan(multiprocess) # choose an appropriate plan from future package
    #plan(multicore)
    n_workers = parallel::detectCores() - 2
    plan(multisession, workers = n_workers )
    cat( "Running simulation in parallel with", n_workers, "workers\n" )

    factors$res <- future_pmap(factors, .f = file_saving_sim,
                          .options = furrr_options(seed = NULL),
                          .progress = TRUE )

} else {
    # Run sequentially
    factors$res <- pmap(factors, .f = file_saving_sim )
}

tictoc::toc()


##### Clean up simulation results for saving #####

factors
head( factors$res[[1]] )

sim_results <-
    factors %>%
    unnest(cols = res)


sim_results

# Cut apart the results and error messages
sim_results$sr = rep( c("res","err"), nrow(sim_results)/2)
sim_results = pivot_wider( sim_results, names_from = sr, values_from = res )
sim_results

# Tally errors
a = sim_results$err
lgc = map_lgl( a, is.null )
table( lgc )
if ( all( lgc ) ) {
    cat( "No errors detected\n" )
} else {
    cat( "Errors detected\n" )
    lgc = map_dbl( a, length )
    print( table( lgc ) )
}

saveRDS( sim_results, file=paste0( result_dir, "simulation_results.rds" ) )

cat( "Script finished\n" )


