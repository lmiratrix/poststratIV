

# This file has code to generate fake data for simulation studies.
#
# The DGP allows for mixing compliers and nevertakers, etc.

library( tidyverse )

#### Data making process ####

if ( FALSE) {

    # parameters
    N = 1000
    shift = -3
    tau = 3

}






#' Randomize into tx and control and observe outcomes
#'
#' Useful for running simulation studies.
#'
#' @param df Dataframe of potential outcomes.  In particular has
#'   columns Y0, Y1, and S0, S1.
#' @param p Proportion of units that are treated.
#'
#' @export
rand.exp = function( df, p=0.3 ) {
    df = dplyr::mutate( df,
                 Z = as.numeric( sample( dplyr::n() ) <= dplyr::n()*p ),
                 S = ifelse( Z, S1, S0 ),
                 Yobs = ifelse( Z, Y1, Y0 ) )
    df
}




#' Make stratified noncompliance dataset
#'
#' Make dataset with different strata with different numbers of
#' compliers in each strata.
#'
#' @param pstrat Proportion of units in each strata.
#' @param pi_c Proportion of compliers in each strata.
#' @param pi_n Proportion of never-takers in each stata, defaults to
#'   remainder making no always-takers.
#' @param Ybar0 baseline mean outcome of each strata.
#' @param sd0 Within-group individual residual variance.
#' @param tau Complier average causal effect for each strata.
#' @param at_shift Always-takers will have their potential outcomes
#'   shifted by this much, relative to compliers.
#' @param nt_shift Never-takers will have their potential outcomes
#'   shifted by this much, relative to compliers.
#' @return Dataframe of individual level data.  X denotes strata
#'   membership (categorical).
#'
#' @export
make.dat.tuned = function( N,
                           pstrat = rep( 0.25, 4 ),
                           pi_c = c( 0, 0.01, 0.12, 0.30 ),
                           pi_n = 1 - pi_c,
                           Ybar0 = c( 2, 5, 7, 10 ),
                           sd0 = 3,
                           tau = c( 0, 1, 2, 3 ),
                           at_shift = -1,
                           nt_shift = 1,
                           exclusion_impact = 0,
                           maintain_tau = TRUE ) {

    K = length( pstrat )

    stopifnot( length(pi_c) == K )
    stopifnot( length(pi_n) == K )
    stopifnot( length(Ybar0) == K )
    stopifnot( length(tau) == K )
    stopifnot( length(pstrat) == K )
    stopifnot( all( pi_c + pi_n <= 1 ) )
    stopifnot( all( pi_c >= 0 ) )
    stopifnot( all( pi_n >= 0 ) )

    if ( length( sd0 ) == 1 ) {
        sd0 = rep( sd0, K )
    }

    # Covariate
    X = sample( 1:K, N, replace=TRUE, prob=pstrat )

    # Make prognostic covariate
    # X2 = sample( 1:K2, N, replace=TRUE )

    # Strata membership
    strat = runif(N)
    # ordered as at, comp, nt

    # Take tx under control?
    pi_a = 1 - pi_c - pi_n
    S0 = 0 + (strat <= pi_a[X] )

    # Take tx under tx?
    S1 = 0 + (strat <= pi_a[X] + pi_c[X] )

    #Y0 = beta*X2 + Ybar0[X] + (1-S0)*(1-S1)*nt_shift + S0*S1*at_shift + rnorm( N, sd = sd0[X] )
    Y0 = Ybar0[X] + (1-S0)*(1-S1)*nt_shift + S0*S1*at_shift + rnorm( N, sd = sd0[X] )

    # If complier, Y1 = Y0 + tau, otherwise Y1 = Y0
    Y1 = Y0 + ifelse( !S0 & S1, tau[X], exclusion_impact )

    df = tibble::tibble( X = X, Y0 = Y0, Y1 = Y1, S0 = S0, S1 = S1,
                 complier = 0 + ( S0 == 0 & S1 == 1 ))

    df$X = paste0( "X", df$X )


    return( df )
}





#' Describe a dataset with compliance
#'
#' Given a dataframe that has both potential outcomes (e.g., a simulated
#' complete dataset), calculate some summary measures of interest.
#'
#' @return Dataframe with several summary measures.  sd0, the variance of the Y0
#'   in the dataframe.
#'
#' @export
describe_sim_data = function( df ) {
    #dfG = df %>%
    #    group_by( X ) %>%
    #    nest()
    #dfG

    #gsum = set_names( dfG$data, nm=dfG$X ) %>%
    #    map_df( IV.est.oracle, include_est=TRUE, .id="X" ) %>%
    #    arrange( X )
    #gsum
    gsum = df %>% group_by( X ) %>%
        summarise( n = n(),
                   n_comp = sum( complier ),
                   Ybar0 = mean( Y0 ),
                   ITT = mean( Y1 - Y0 ),
                   pi = mean( complier ),
                   nt = mean( S0 == 0 & S1 == 0 ),
                   at = mean( S0 == 1 & S1 == 1 ),
                   LATE = ITT / pi,
                   sd0 = sd( Y0 ) )


    # Aggregate
    gagg = df %>%
        summarise( n = n(),
                   n_comp = sum( complier ),
                   Ybar0 = mean( Y0 ),
                   ITT = mean( Y1 - Y0 ),
                   pi = mean( complier ),
                   nt = mean( S0 == 0 & S1 == 0 ),
                   at = mean( S0 == 1 & S1 == 1 ),
                   LATE = ITT / pi,
                   sd0 = sd( Y0 ) )
    gagg$X = "ALL"

    gsum = dplyr::bind_rows( gsum, gagg )

    gsum = dplyr::mutate( gsum, per = 200*n / sum(n),
                   per_c = 200*n_comp / sum(n_comp) )

    gsum

}



#' Summarize simulated dataset
#'
#' Calculate how predictive the covariate X is in terms of outcome and of
#' compliance.
#'
#' @param dd Dataframe of simulation data with three columns: Yobs,
#'   complier (a flag of being a complier) and X (a covariate).
#'
#' @return Dataframe with R2_y and R2_comp from two regressions.
#'
#' @export
summarize_sim_data = function( dd ) {
    M0 = lm( Yobs ~ X, subset( dd, Z == 0 ) )
    sm <- summary( M0 )
    sm$r.squared

    M2 = lm( complier ~ X, subset( dd, Z == 0 ) )
    sm2 <- summary( M2 )
    sm2$r.squared

    tibble( R2_y = sm$r.squared, R2_comp = sm2$r.squared )
}


##### Driver sim data method ######



#' Simulate data for non-compliance in an RCT
#'
#' Make data given set parameters.  Idea is to have standardized
#' outcomes regardless of parameters, mostly.
#'
#' @param nt_shift shifts the never takers mean outcome so they are
#'   different from compliers.
#' @param sd0 baseline Y0 standard deviation not accounting for shift
#'   of NT and AT and treatment impacts.
#' @param pi_c Overall proportion of compliers
#' @param pi The relative proportion of compliers in each strata (this
#'   will be rescaled to achieve overall pi_c proportion compliers).
#' @param one_sided If FALSE, split non-compliers into AT and NT.
#' @param frac_nt The fraction of non-compliers that are never takers
#'   (can vary by strata).
#' @param pred_comp = FALSE means all strata have same proportion of
#'   compliers
#' @param pred_Y = FALSE means all strata have same Y0 means (for
#'   compliers).
#' @param perfect_X = TRUE means make X4 all compliers, and X1, X2,
#'   and X3 random mix of the noncompliers.
#' @param scaled_C If a number, pack compliers into the strata with
#'   proportions equal to a geometric series with, for X4, X3, X2 and
#'   X1 relative proportion 1, scaled_C, scaled_C^2, scaled_C^3 such
#'   that the sum equals the target desired pi_c.
#' @param exclusion_impact Treatment impact for noncompliers (i.e.
#'   violation of the exclusion restriction).
#' @param pstrat The relative sizes of strata from X1 to X4.
#'   Overridden if params is not null
#' @param params Dataframe with pstrat, pi, Ybar0, and tau as default
#'   values for the strata.  Will overwrite the passed values of
#'   these, if params is not null.
#' @param include_POs Should the potential outcomes be included in the
#'   final returned data?  If FALSE, then you get data that would look
#'   like what you would get in the field.  TRUE means you have the
#'   underlying truth to calculate actual (oracle) effects.
#'
#' @export
make_dat = function( N, pi_c = NULL,
                     frac_nt = 1/3,
                     one_sided = TRUE, at_shift = 0, nt_shift = 0, sd0 = 1,
                     pred_comp = "yes", pred_Y = "yes", het_tx = "yes",
                     exclusion_impact = 0,
                     perfect_X = FALSE,
                     scaled_C = NULL,
                     verbose = FALSE,
                     pstrat = rep( 0.25, 4 ),
                     pi = c( 0.01, 0.05, 0.55, 0.95 ),
                     Ybar0 = c( 4, 3, 2, 1 ),
                     tau = c( 0, 1, 2, 3 ),
                     params = NULL,
                     include_POs = TRUE ) {

    # Copy over params, overriding the other passed parameters
    if ( !is.null( params ) ) {
        pstrat = params$pstrat
        pi = params$pi
        Ybar0 = params$Ybar0
        tau = params$tau
    }

    if ( is.null( pi_c ) ) {
        pi_c = weighted.mean( pi, pstrat )
    }

    K = length( pstrat )
    stopifnot( length( Ybar0 ) == K )
    stopifnot( length( pi ) == K )
    stopifnot( length( tau ) == K )

    # These combos of parameters don't make sense.
    stopifnot( perfect_X == FALSE || is.null(scaled_C) )
    stopifnot( pred_comp == "yes" || is.null(scaled_C) )


    # The strata-level treatment impacts
    if ( het_tx == "no" ) {
        tau = rep( weighted.mean( tau, pstrat * pi_c), K )
    }

    if ( pred_comp == "no" ) {
        pi_c = rep( pi_c, K )
    } else if ( !is.null( scaled_C ) ) {
        scl = c( scaled_C^3, scaled_C^2, scaled_C, 1 )
        wt = pstrat * scl
        p = pi_c / sum(wt)
        stopifnot( p <= 1 )
        pi_c = p * scl
    } else {
        # rescale complier proportions to hit overall target.
        pi_cur = weighted.mean( pi, pstrat )
        pi_c = pi * (pi_c / pi_cur)
    }

    stopifnot( all( pi_c <= 1 ) )
    stopifnot( all( pi_c >= 0 ) )

    Ybar0 = Ybar0
    grand_mean = weighted.mean( Ybar0, w=pstrat )
    if ( pred_Y == "no" ) {
        Ybar0 = rep( grand_mean, length( Ybar0 ) )
    } else {
        # deflate residual variation to account for group
        # variation
        sd0 = sd0^2 - weighted.mean( (Ybar0 - grand_mean)^2,
                                     w = pstrat )
        if ( sd0 <= 0 ) {
            browser()
        }
        stopifnot( sd0 >= 0 )
        sd0 = sqrt( sd0 )
    }

    # Print some stuff about our DGP, if desired.
    if ( verbose ) {
        scat( "\nGen data (comp=%s Y=%s het=%s): N = %d, sd0 = %.3f\n",
              pred_comp, pred_Y, het_tx, N, sd0 )
        print( tibble( pstrat = pstrat,
                       pi_c = pi_c,
                       Ybar0 = Ybar0,
                       nt_shift = nt_shift,
                       tau = tau
        ))
    }

    # Note: if not one-sided, create space for always takers
    if ( one_sided ) {
        pi_nt = 1 - pi_c
    } else {
        #browser()
        #cat( "two-sided\n" )
        pi_nt = frac_nt * (1-pi_c)
    }

    df = make.dat.tuned( N = N,
                         pstrat = pstrat,
                         pi_c = pi_c,
                         pi_n = pi_nt,
                         Ybar0 = Ybar0,
                         sd0 = sd0,
                         nt_shift = nt_shift,
                         at_shift = at_shift,
                         exclusion_impact = exclusion_impact,
                         tau = tau )

    # Make a strata that is all compliers, and the other strata almost
    # entirely not compliers.
    if ( perfect_X ) {
        df$X = sample( 1:3, nrow(df), replace=TRUE )
        df$X[df$complier == 1] = 4
        df$X = paste0( "X", df$X )

        # Move 1 complier to each of the other strata
        comps = which( df$complier == 1 )
        df$X[comps[1:3]] = c( "X1", "X2", "X3")
    }

    df = rand.exp( df )

    if ( !include_POs ) {
        df = df |> dplyr::select( -Y0, -Y1, -S1, -S0 )
    }
    return( df )
}






