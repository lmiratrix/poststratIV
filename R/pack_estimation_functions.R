
# This script has functions to implement different versions of IV estimation
# (both stratified and not).

library( "AER" )


#### Utility functions ####

#' For pretty printing
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#### The estimators (including oracle) ####



#' Parameter calculation and oracle estimation (For simulation studies)
#'
#' This function does two things:
#'
#' (1) Calculate some simulation statistics (sample parameter values) This is
#' done with the full schedule of potential outcomes, something we would not see
#' empirically.
#'
#' (2) (If include_est = TRUE) estimate (using an oracle estimator) the CACE.
#' The oracle estimator is the difference in means on the compilers only.
#'
#' (this method is for simulation code).
#'
#' @param data Dataframe with named Y0, Y1, S0, S1 and Z columns.
#'
#' @export
#' @importFrom dplyr n
IV.est.oracle = function( data, include_est=FALSE  ) {

    # First get true parameters for our data for our different princ strata (we
    # are a simulation so we can!)
    stats = data |> dplyr::group_by( S0, S1 ) |>
        dplyr::summarize( Y.bar0 = mean( Y0 ),
                          Y.bar1 = mean( Y1 ),
                          p = n() / nrow(data) )
    comp = dplyr::filter( stats, S0 == 0 & S1==1 )

    # True finite sample ITT and LATE and proportion compliers
    if ( nrow( comp ) == 1 ) {
        res = data.frame( ITT = mean(data$Y1 - data$Y0),
                          LATE = comp$Y.bar1 - comp$Y.bar0,
                          pi = comp$p )
    } else {
        # No compliers, so fail
        res = data.frame( ITT = 0, LATE = NA, pi = 0 )
    }
    res$n = nrow(data)

    # If desired, run an oracle estimator of a difference in means on just the
    # compliers (if we somehow knew who they were)
    if ( include_est ) {
        ss = dplyr::filter( data, S0 == 0, S1 == 1 )
        stats = ss |> dplyr::group_by( Z ) |> dplyr::summarize( Y.bar = mean( Yobs ),
                                                                var.Y = var( Yobs ),
                                                                n = n() )
        if ( nrow(stats) > 1 ) {
            res = mutate( res,
                          LATE.oracle = diff( stats$Y.bar ),
                          SE.oracle = with(stats, sqrt( sum( var.Y / n ) ) ) )
        } else {
            res = mutate( res,
                          LATE.oracle = NA,
                          SE.oracle = NA )
        }
    }

    res
}


#' Implement different versions of classic IV estimation.
#'
#' @param Yobs, S, Z : column names of outcome, treatment receipt and
#'   treatment assignment (instrument), respectively.
#'
#' @return One row tibble of results.  If not both treatment and
#'   control units present, return tibble with NAs for estimated
#'   impacts, etc.
#' @export
IV.est = function( data, Yobs="Yobs", S="S", Z="Z",
                   tolerance = 1/10000 ) {

    stopifnot( "S" != Yobs && "S" != Z )
    stopifnot( "Yobs" != S && "Yobs" != Z )
    stopifnot( "Z" != Yobs && "Z" != S )

    data$Yobs = data[[Yobs]]
    data$S = data[[S]]
    data$Z = data[[Z]]

    # Use observed outcome data to estimate our parameters
    # First get summary statistics
    stats = data |> dplyr::group_by( Z ) |>
        dplyr::summarize( Y.bar = mean( Yobs ),
                          p.S = mean(S),
                          var.Y = var( Yobs ),
                          var.S = var( S ),
                          cov.Y.s = cov( Yobs, S ),
                          n = n() )

    # Need tx and co groups
    if ( nrow(stats) != 2 ) {
        return( tibble( ITT.hat = NA,
                        pi.hat = NA,
                        LATE.hat = NA,
                        SE.ITT = NA,
                        SE.pi = NA,
                        cov.ITT.pi = NA,
                        SE.wald = NA,
                        SE.delta = NA,
                        n = nrow(data)
        ))
    }
    #stopifnot( nrow( stats ) == 2 )

    pi.hat = diff( stats$p.S )
    if ( abs(pi.hat) < tolerance ) {
        #if ( abs(pi.hat) > 0.000000000001 ) {
        #    browser()
        #}
        pi.hat = 0
    }

    ITT.hat = diff( stats$Y.bar )

    # Store estimates
    res = data.frame( ITT.hat=ITT.hat,
                      pi.hat = pi.hat,
                      LATE.hat = ITT.hat / pi.hat )

    # Calc primary SEs of the direct estimands
    res = mutate( res,
                  SE.ITT = with( stats, sqrt( sum( var.Y / n ) ) ),
                  SE.pi = with( stats, sqrt( sum( var.S / n ) ) ),
                  cov.ITT.pi = with( stats,  sum( cov.Y.s / n )  ) )

    # Calculate SEs, including delta method SE (dropping correlation term)
    res = mutate( res,
                  SE.wald = SE.ITT / pi.hat,
                  SE.delta = sqrt( SE.ITT^2 / pi.hat^2 + ITT.hat^2 * SE.pi^2 / pi.hat^4 - 2*ITT.hat * cov.ITT.pi / pi.hat^3 )
    )

    res$n = nrow( data )

    # clean up if pi.hat = 0
    if ( res$pi.hat == 0 ) {
        res$LATE.hat = NA
        res$SE.wald = NA
        res$SE.delta = NA
    }

    # and done!
    res
}





#### Stratified estimator ####



#' Aggregate strata using specified weighting and calculate associated
#' weighted-average Standard Error estimates as well.
#'
#' @param gsum Data table of strata-level estimates (point estimates
#'   and standard errors, etc)
#' @param Weighting: normal means weight by estimated number of
#'   compliers, double means square of that, precision means 1/SE^2.
#' @return Aggregated result
#'
#' @export
aggregate_strata <- function( gsum, weighting = c( "normal", "double", "precision" ) ) {


    if ( nrow( gsum ) == 0 ) {
        return( tibble( Xblk = "IV_w",
                        LATE.hat = NA,
                        SE.wald = NA,
                        SE.delta = NA,
                        n = 0,
                        n_comp = 0,
                        pi.hat = NA,
                        ITT.hat = NA ) )
    }

    weighting = match.arg(weighting)

    if ( weighting == "double" ) {
        gsum$weight = gsum$n_comp^2
    } else if ( weighting == "normal" ) {
        gsum$weight = gsum$n_comp
    } else if ( weighting == "precision" ) {
        gsum$weight = 1 / gsum$SE.wald^2
    } else {
        stop( "Invalid weighting for strata aggregation" )
    }

    # Normalizing constant
    W_c = sum(gsum$weight)^2

    IV_w = gsum %>% summarise( Xblk = "IV_w",
                               LATE.hat = weighted.mean( LATE.hat, w=weight ),
                               SE.wald = sqrt( sum( (weight^2 / W_c) * SE.wald^2 ) ),
                               SE.delta = sqrt( sum( (weight^2 / W_c) * SE.delta^2 ) ),
                               total_n = sum(n),
                               n_comp = sum(n_comp),
                               pi.hat = weighted.mean( pi.hat, w=n ),
                               ITT.hat = weighted.mean( ITT.hat, w=n ),
                               k = nrow(gsum) ) %>%
        dplyr::rename( n = total_n )

    if ( is.nan( IV_w$LATE.hat ) ) {
        IV_w$LATE.hat = NA
        IV_w$SE.wald = NA
        IV_w$SE.delta = NA
    }

    IV_w
} # end aggegate_strata



get_anova_F <- function( data ) {
    s_anova = purrr::quietly( anova )
    res <- s_anova( lm( S ~ Z, data=data ) )

    res$result$`F value`[[1]]
}



# Calculate the IVa estimator using the strata-specific statistics.
calc_IVa <- function( gsum ) {
    W = sum(gsum$n)^2
    gsum$.weight = gsum$n

    gsum2 <- gsum %>%
        summarise( ITT.hat.ps = weighted.mean( ITT.hat, w=.weight ),
                   pi.hat.ps = weighted.mean( pi.hat, w=.weight ),
                   SE.ITT.ps = sqrt(sum(SE.ITT^2*n^2)/W),
                   SE.pi.ps = sqrt(sum(SE.pi^2*n^2)/W),
                   cov.ITT.pi.ps = sum(cov.ITT.pi*n^2)/W,
                   n = sum(n),
                   n_comp = sum(n_comp) )

    stopifnot( nrow(gsum2)==1 )
    if ( is.na( gsum2$ITT.hat.ps ) ) {
        stop()
    }

    IV_a = gsum2 |> mutate( Xblk = "IV_a",
                            LATE.hat = ITT.hat.ps/pi.hat.ps,
                            SE.wald = SE.ITT.ps / abs(pi.hat.ps),
                            SE.delta = sqrt( SE.ITT.ps^2 / pi.hat.ps^2 +
                                                 (ITT.hat.ps/pi.hat.ps)^2 * SE.pi.ps^2 / pi.hat.ps^2 -
                                                 2 * (ITT.hat.ps/pi.hat.ps) * cov.ITT.pi.ps / pi.hat.ps^2 ),
                            n = n,
                            n_comp = n_comp,
                            pi.hat = pi.hat.ps,
                            ITT.hat = ITT.hat.ps,
                            SE.ITT = SE.ITT.ps,
                            SE.pi = SE.pi.ps,
                            cov.ITT.pi=cov.ITT.pi.ps)

    if ( is.infinite( IV_a$LATE.hat ) ) {
        IV_a$LATE.hat = NA
        IV_a$SE.wald = NA
        IV_a$SE.delta = NA
    }

    IV_a
}




#' This will group units based on X and then estimate the LATE within
#' each subgroup. It then averages the subgroups with weight
#' proportional to the estimated proportion of compliers in each
#' group.
#'
#' STRAT_cw: Estimate LATE in each strata, then take the weighted
#' mean, typically weighting by estimated number of compliers.
#'
#' STRAT2: Estimate the post-stratified ITT and pi, take the ratio to
#' get overall LATE estimate.
#'
#' @param Yobs, S, Z : column names of outcome, treatment receipt and
#'   treatment assignment (instrument), respectively.
#' @param strat_var Name of column in data that holds a categorical
#'   variable that we should stratify on.
#' @param drop_without_warning If TRUE will drop blocks that are
#'   inestimable due to not having at least 2 tx and 2 co units.  If
#'   FALSE will drop, but also throw warning.
#' @param tidy_table TRUE means only return subset of results in the
#'   table (primary results of interest)
#' @export
#' @importFrom rlang set_names
#' @importFrom dplyr group_by filter mutate arrange summarise
#'   left_join bind_rows relocate n
#' @importFrom dplyr mutate
#' @importFrom rlang ensym
IV.est.strat = function( data,  Yobs="Yobs", S="S", Z="Z", strat_var ="X",
                         DSS_cutoff = 0.02,
                         drop_without_warning = FALSE,
                         include_blocks = TRUE,
                         include_FSS = TRUE,
                         return_strata_only = FALSE,
                         tidy_table = FALSE ) {


    stopifnot( "S" != Yobs && "S" != Z )
    stopifnot( "Yobs" != S && "Yobs" != Z )
    stopifnot( "Z" != Yobs && "Z" != S )

    data$Yobs = data[[Yobs]]
    data$S = data[[S]]
    data$Z = data[[Z]]

    # give canonical name to our passed variable
    data$Xblk = as.character( data[[strat_var]] )
    stopifnot( !is.null(data$Xblk) )
    if ( any( is.na( data$Xblk ) ) ) {
        if ( "NA" %in% data$Xblk ) {
            stop( glue::glue( "Stratification covariate has missing values. Automatically converting NAs to 'NA' in column '{strat_var}' collides with existing 'NA' in blocking covariate." ) )
        } else {
            warning( "Stratification covariate has missing values. These will be formed into a separate strata called 'NA', but we recommend doing this manually as a pre-processing step." )
        }
        data$Xblk[ is.na( data$Xblk ) ] = "NA"
    }

    dataG = data |>
        dplyr::group_by( Xblk ) |>
        tidyr::nest()

    gsum = rlang::set_names( dataG$data, nm=dataG$Xblk ) |>
        purrr::map_df( IV.est, .id="Xblk" )



    # Zero out empty estimated strata
    empty = gsum$pi.hat <= 0 | !is.finite( gsum$LATE.hat )
    gsum$LATE.hat[ empty  ] = 0
    gsum$SE.wald[ empty ] = 0
    gsum$SE.delta[ empty ] = 0

    # Estimate the number of compliers in each strata
    gsum = mutate( gsum,
                   n_comp = pi.hat*n ) |>
        arrange( Xblk )

    # Drop all strata without at least 2 Tx and Co units
    nr = nrow(gsum)
    gsum = dplyr::filter( gsum, !is.na( SE.ITT ) )
    nr_post = nrow(gsum)
    if ( !drop_without_warning && nr > nr_post ) {
        warning( glue::glue( "{nr - nr_post} blocks dropped (out of {nr}) due to fewer than 2 Tx or 2 Co units" ) )
    }

    if ( return_strata_only ) {
        return( gsum )
    }


    #### The estimators ####


    # IV_w: The weighted average of the strata-level estimates
    IV_w <- gsum |>
        dplyr::filter( pi.hat != 0 ) |>
        aggregate_strata()
    #if ( IV_w$SE.delta > 10 ) {
    #    browser()
    #}

    # The "drop 0 or less" estimator
    res_DSS_0 <- gsum |>
        dplyr::filter( pi.hat > 0 ) |>
        aggregate_strata()
    res_DSS_0$Xblk = "DSS0"

    # The "drop small strata" estimator
    res_DSS = gsum |>
        dplyr::filter( pi.hat > DSS_cutoff ) |>
        aggregate_strata()
    res_DSS$Xblk = "DSS"


    # IV_a: The other version of stratification: First calculate
    # post-stratified ITT.hat, SE.ITT, pi.hat, and SE.pi, and then
    # calculate LATE with the post-stratified estimates of ITT and pi.

   IV_a <- calc_IVa( gsum )


    # "Drop Small F", or "Testing-F" estimator
    res_DSF = NULL
    if ( include_FSS ) {
        Fvals <- data |>
            dplyr::nest_by( Xblk ) |>
            summarise( Fval = get_anova_F( data ),
                       .groups = "drop" )
        res_DSF <- gsum |>
            left_join( Fvals, by="Xblk" ) |>
            dplyr::filter( Fval > 10 )
        if ( nrow( res_DSF ) == 0 ) {
            warning( "All strata dropped for DSF method" )
            res_DSF = tibble( LATE.hat = NA )
        } else {
            res_DSF <- aggregate_strata( res_DSF )
        }
        res_DSF$Xblk = "DSF"
    }


    # The "double weight" estimator
    # res_DDW = gsum |>
    #     dplyr::filter( pi.hat > 0 ) |>
    #     aggregate_strata( weighting = "double" )
    # res_DDW$Xblk = "DDW"

    # The "double weight" (precision-weighted, actually) estimator
    res_DDW = gsum |>
        dplyr::filter( pi.hat > 0 ) |>
        aggregate_strata( weighting = "precision" )
    res_DDW$Xblk = "PWIV"

    # No stratification as baseline
    res_unstrat = IV.est( data )
    res_unstrat = mutate( res_unstrat,
                          Xblk = "UNSTRAT",
                          n_comp = n * pi.hat )

    gsum$k = 1
    res_unstrat$k = nrow(gsum)
    IV_a$k = nrow(gsum)

    if ( !include_blocks ) {
        gsum = NULL
    }

    res = bind_rows( gsum, res_unstrat, IV_w, res_DSS_0, IV_a, res_DSS, res_DDW, res_DSF ) |>
        relocate( Xblk )

    res$empty = sum( empty )

    if ( any( is.nan( res$LATE.hat ) ) || any( is.infinite( res$LATE.hat ) ) ) {
        stop( "Uneexpected NaN or Infinite value" )
    }

    if ( tidy_table ) {
        res = tidy_IV_table( res )
    }

    res
}




#' Estimate stratifying on a continuous variable X
#'
#' This just cuts provided continuous covariate X into 4 bins and does the post-stratification on the
#' resulting categorical variable.
#'
#' @param df Dataframe with standard outcome columns.
#' @param strat_var Name of column in data that holds a continuous covariate.
#'
#' @export
IV.est.strat.cont = function( df, strat_var="X" ) {
    df$Xblk = as.character( cut( df$X, 4 ) )
    IV.est.strat( df, "Xblk" )
    #    IV.est( df )
}


#' Tidy up IV results table
#'
#' This simply selects some more relevant columns from the extensive
#' output of the main method.
#'
#' @noRd
tidy_IV_table = function( res ) {
    dplyr::select( res, Xblk, ITT.hat:LATE.hat,
                   SE.wald, k, n, n_comp, empty ) |>
        mutate( n_comp = round( n_comp ) )
}
