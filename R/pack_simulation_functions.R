
# Simulation function library



# This simulation study focuses on the LATE given a DGP that has CACE vary by
# proportion of complier.
#
# It examines two things:
# (1) The performance of the point estimate
# (2) The performance of the standard error estimate
#     - do we get good coverage?
#     - Is the estimate of the SE estimate more or less unstable?






# These are estimated from the baseline empirical data
# default_params = tibble::tribble(
#     ~pstrat,  ~pi, ~Ybar0,  ~tau, ~sd0,
#     0.07, 0.09,   0.08,  0.08, 0.28,
#     0.54, 0.03,   0.11, -0.41, 0.31,
#     0.29, 0.17,   0.18, -0.03, 0.38,
#     0.09, 0.02,   0.26, -0.15, 0.44
# )

# THese are an earlier sim round
# params = tibble::tribble(
#     ~pstrat,  ~pi, ~Ybar0,  ~tau,
#     0.40, 0.025,   0.00,    0.0,
#     0.25, 0.05,    0.25,    0.1,
#     0.20, 0.10,    0.50,    0.2,
#     0.15, 0.80,    0.75,    0.4
# )


#' Default simulation parameters for simulation code
#'
#' @return Small tibble with simulation parameters that can be passed
#'   to simulation functions.
#'
#' @export
default_sim_params <- function( ) {


    sim_params = tibble::tribble(
        ~pstrat,  ~pi, ~Ybar0,   ~tau,
        0.40,   0.01,    4,    0.0,
        0.45,   0.05,    3,    0.5,
        0.10,   0.55,    2,    1.0,
        0.05,   0.95,    1,    1.5
    )

    return( sim_params )
}

#
# params = tibble::tribble(
#     ~pstrat,  ~pi, ~Ybar0,   ~tau,
#     0.25,   0.05,    4,    0.0,
#     0.25,   0.10,    3,    0.1,
#     0.25,   0.45,    2,    0.2,
#     0.25,   0.90,    1,    0.4
# )






#### Testing data generator function #####

# Print out statistics about the data generated via make_dat

if ( FALSE ) {

    #    params$tau = 10 * params$tau

    grd = expand_grid( pred_comp = c("no","yes"),
                       pred_Y = c("no", "yes" ),
                       het_tx = c("no", "yes" ) )

    grd$data = pmap( grd, make_dat, N = 100000, nt_shift = 1,
                     sd0 = 1, pi_c = 0.10, params = default_sim_params() )

    grd = tidyr::unnest( grd, cols=data )
    grd
    table( S0=grd$S0, S1=grd$S1 ) / nrow(grd)

    sgrd <- grd |> group_by( pred_comp, pred_Y, het_tx ) |>
        summarise( n = n(),
                   sdY0 = sd( Y0 ),
                   sdY1 = sd( Y1 ),
                   pcomp = mean( complier ),
                   ptx = mean( Z ) )

    sgrd = pivot_longer( sgrd, cols = sdY0:ptx, names_to="var", values_to="value"  )

    head( sgrd )
    sgrd$value = round( 2*sgrd$value, digits = 2 )/2
    ggplot( sgrd, aes( pred_comp, value, col=pred_Y, pch=het_tx ) ) +
        facet_wrap( ~ var, scales="free", nrow=1 ) +
        geom_jitter( width = 0.3, height=0)


    datA = make_dat( 10000, pi_c = 0.10, pred_comp = "no", het_tx = "yes", sd0 = 1, nt_shift = -1,
                     params = sim_params)
    describe_sim_data(datA)

    datB = make_dat( 10000, pi_c = 0.10, pred_comp = "yes", het_tx = "yes", sd0 = 1, nt_shift = -1,
                     params = sim_params )
    describe_sim_data(datB)
}



##### 2SLS Baseline ######

est.2SLS = function( df ) {
    # 2SLS
    reg_iv0 <- AER::ivreg( Yobs ~ X + S|X + Z, data = df)
    cc = coef( reg_iv0 )[["S"]]
    sm = summary( reg_iv0 )
    SE = NA
    if ( "S" %in% rownames( sm$coefficients ) ) {
        SE = sm$coefficients["S", "Std. Error"]
    }
    tt = tibble( Xblk = "2SLS",
                 LATE.hat = cc,
                 SE.wald = SE,
                 SE.delta = SE,
                 n = nrow(df) )

    # reg_null <- ivreg( Yobs ~ S|Z, data = df)
    # cc = coef( reg_null )[["S"]]
    # sm = summary( reg_null )
    # SE = sm$coefficients["S", "Std. Error"]
    # tt2 = tibble( Xblk = "2SLS_unadj",
    #              LATE.hat = cc,
    #              SE.wald = SE,
    #              SE.delta = SE,
    #              n = nrow(df) )

    tt
}



##### The simulation functions #####



one_run = function( N, nt_shift, sd0, pi_c,
                    pred_comp, pred_Y, het_tx,
                    perfect_X = FALSE, scaled_C = NULL,
                    shuffle = FALSE,
                    one_sided = TRUE,
                    params = default_sim_params(),
                    data_only = FALSE, ...  ) {

    stopifnot( pred_comp %in% c( "no", "yes" ) )
    stopifnot( pred_Y %in% c( "no", "yes" ) )
    stopifnot( het_tx %in% c( "no", "yes" ) )

    df = make_dat( N = N, nt_shift = nt_shift, sd0 = sd0, pi_c = pi_c,
                   pred_comp = pred_comp, pred_Y = pred_Y, het_tx = het_tx,
                   one_sided = one_sided,
                   perfect_X = perfect_X, scaled_C = scaled_C, params = params, ... )

    if ( shuffle ) {
        df$X = sample( df$X )
    }

    if ( data_only ) {
        return( df )
    }

    IVest = IV.est.strat(df)

    g = IVest$Xblk
    IVest$pi_zero <- 0
    if ( any( is.na( IVest$pi.hat ) ) ) {
        IVest$pi_zero = 1 # Undefined in strata means pi_zero is de-facto 0?
    } else if( any( IVest$pi.hat[startsWith(IVest$Xblk, "X")] == 0) ) {
        IVest$pi_zero <- 1
    }

    IVest = dplyr::filter( IVest, !( Xblk %in% c( "X1", "X2", "X3", "X4" ) ) )

    orc = IV.est.oracle( df, include_est = TRUE )
    orc = dplyr::rename( orc,
                  LATE.hat = LATE.oracle,
                  SE.wald = SE.oracle,
                  pi.hat = pi )
    orc$Xblk = "Oracle"
    orc$n_comp = orc$n * orc$pi
    orc$pi_zero<-0


    tt = est.2SLS( df )


    IVest = bind_rows( IVest, orc, tt )

    IVest$ITT = orc$ITT
    IVest$LATE = orc$LATE
    IVest$pi = orc$pi
    IVest = rename( IVest,
                    method = Xblk )
    IVest
}



#' Run a simulation of noncompliance in an RCT
#'
#' @param reps Number of replications
#' @param N Number of observations
#' @param nt_shift Shift in the mean of the outcome for the
#'   never-taker group.
#' @param sd0 Standard deviation of the outcome
#' @param pi_c Proportion of compliers
#' @param pred_comp Is the covariate predictive of compliance?
#' @param pred_Y Is the covariate predictive of the outcome?
#' @param het_tx Is the treatment effect heterogeneous?
#' @param one_sided Is compliance one-sided (TRUE) or two-sided
#'   (FALSE)?
#' @param scaled_C passed to DGP function
#' @param perfect_X passed to DGP function
#' @param shuffle Should the X variable be shuffled to break all
#'   connection to anything else?
#' @param params Simulation parameters as a tibble (see, e.g.,
#'   default_sim_params()).
#' @param seed Seed for the random number generator
#' @param ... Additional arguments passed to the DGP function
#'
#' @return A tibble with the results of the simulation.
#'
#' @export
run_sim <- function( reps,
                     N, nt_shift, sd0, pi_c,
                     pred_comp, pred_Y, het_tx,
                     one_sided = TRUE,
                     scaled_C = NULL,
                     perfect_X = FALSE,
                     shuffle = FALSE,
                     params = default_sim_params(),
                     seed = NULL, ... ) {

    options(list(dplyr.summarise.inform = FALSE))

    if (!is.null(seed)) set.seed(seed)

    runs <-
        purrr::map( 1:reps,
                    function(x) {
                        one_run( N=N, nt_shift = nt_shift,
                                 sd0 = sd0, pi_c = pi_c,
                                 one_sided = one_sided,
                                 pred_comp = pred_comp,
                                 pred_Y = pred_Y,
                                 het_tx = het_tx,
                                 perfect_X = perfect_X,
                                 scaled_C = scaled_C,
                                 shuffle = shuffle,
                                 params = params,
                                 ... )
                    } ) |>
        bind_rows( .id="runID" )

    return( as_tibble( runs ) )
}



