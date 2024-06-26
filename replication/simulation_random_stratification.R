
# This simple simulation shows that if you post-stratify essentially at random
# in the case of low compliance, you gain substantially in your standard errors.

# It indeed appears as if dropping big chunks of folks that we know
# have no ability to help estimate impacts is a boon.  That makes
# sense, but doing it at random?  Yikes!

# Next question is how many chunks is ideal to get the most gains?


library( tidyverse )
library( poststratIV )
library( cli )

cat( "\n\n" )
cli_alert_info( "Running the random stratification simulation\n" )



if ( !exists( "R" ) ) {
    # Number of sim runs
    R = 6
    warning( "Setting R=6 to demo code.  Raise number of simulation replicates for real results" )
    
    # M_CHUNK is a hack to make a bunch of smaller chunks for doing parallel more
    # efficiently.
    M_CHUNK = 2
}


# Where to save results
results_file = here::here( "results/random_strat_sim.rds" )

if ( !dir.exists( dirname( results_file ) ) ) {
    dir.create( dirname( results_file ) )
}




one_run_mini = function( num_strata, nt_shift, one_sided, pi_c = 0.05 ) {

    dat = make_dat( N = 500,
                    one_sided = one_sided, at_shift = -1,
                    pi_c = pi_c,
                    nt_shift = nt_shift,
                    sd0 = 1,
                    pred_comp = "no", pred_Y = "no", het_tx = "no" )

    dat$X = sample( LETTERS[1:num_strata], nrow(dat), replace=TRUE )

    LATE = mean( dat$Y1 - dat$Y0 ) / mean( dat$S1 - dat$S0 )
    LATE

    dat$Y0 = dat$Y1 = NULL
    dat$S0 = dat$S1 = NULL
    dat$complier = NULL

    IVest = IV.est.strat(dat)

    IVest <- IVest %>% filter( Xblk %in% c( "UNSTRAT", "IV_a", "IV_w" ) ) %>%
        rename( method = Xblk )

    IVest$LATE = LATE

    IVest$LATE.hat = round( IVest$LATE.hat, digits=10 )

    if ( num_strata == 1 && (!is.na( IVest$LATE.hat[[1]] ) && IVest$LATE.hat[[1]] != IVest$LATE.hat[[2]] ) ) {
        browser()
    }

    IVest
}



if ( FALSE ) {
    tst = one_run_mini( 10, nt_shift=2, one_sided = FALSE )
    tst


    # Testing 1 strata
    dat = make_dat( N = 500,
                    nt_shift = -1, at_shift = -0.5,
                    pi_c = 0.05, one_sided = FALSE,
                    sd0 = 1,
                    pred_comp = "no", pred_Y = "no", het_tx = "no" )

    head( dat )
    dat$X = "A"
    IVest = IV.est.strat(dat)
    IVest
    table( S0=dat$S0, S1=dat$S1 )

}


sim_mini = function( num_strata, nt_shift, R, one_sided ) {
    rps = rerun( R, one_run_mini( num_strata = num_strata, nt_shift = nt_shift, one_sided = one_sided ) )
    rps = do.call( bind_rows, rps )
    rps
}
#safe_run_sim = safely( sim_mini )


##### Set up simulation #####

dt = expand_grid( num_strata = c( 1, 3, 6, 9, 12 ),
                  nt_shift = c( -1, 0, 1, 2 ),
                  one_sided = c( TRUE, FALSE ),
                  chunk = 1:M_CHUNK)

dt
dt$chunk = NULL


if ( FALSE ) {

    # 1 strat sim check
    dt = filter( dt, num_strata == 1 )
    nrow( dt )
    dt$res = pmap( dt, sim_mini, R = R / M_CHUNK )

}


##### Run the simulation #####

# To run on source, set following to TRUE

if ( TRUE ) {

    parallel::detectCores()

    library(future)
    library(furrr)

    plan(multisession, workers = parallel::detectCores() - 1 )

    # Do the run!  (parallel version)

    tictoc::tic()
    safe_mini = safely( sim_mini )
    res <- future_pmap( dt, safe_mini, R = R / M_CHUNK,
                           .options = furrr_options(seed = NULL),
                           .progress = TRUE )
    res = purrr::transpose( res )
    dt$res = res[[1]]
    dt$err = res[[2]]
    cat( "\n" )
    tictoc::toc()

    # Unpack and save results

    dtT = as_tibble( dt )
    dtT
    dtT = mutate( dtT, no_err = map_lgl(err, is.null) )

    filter( dtT, !no_err )

    dtT = dtT %>% filter( no_err ) %>%
      unnest( cols="res" )
    dtT

    summary( dtT$LATE.hat )

    table( is.na( dtT$LATE.hat ), dtT$method )


    saveRDS( dtT, file = results_file )
}



