
# Analyze the results of the bonus simulation of simulation_pred_C.R
# that looks at an ever-increasingly predictive covariate that only
# predicts compliance status.



library( tidyverse )
library( scales )

rm(list=ls())

theme_set( theme_minimal() )

library( poststratIV )
library( cli )

cat( "\n\n" )
cli_alert_info( "Analyzing the predictive-of-compliance simulation\n" )


source( "plot_helpers.R" )


##### Load and prep simulation data #####

simres1 =  readRDS( here::here( "results_onesided/simulation_results_predC.rds" ) )
simres2 = readRDS( here::here( "results_twosided/simulation_results_predC.rds" ) )
testthat::expect_equal( nrow( simres2 ), nrow( simres1 ) )
simres = bind_rows( one = simres1, two = simres2, .id="sides" )
rm( simres1, simres2 )

simres$err = NULL
simres$seed = NULL
simres$reps = NULL
simres$ID = 1:nrow(simres)
simres <- mutate( simres,
                  one_sided = sides == "one" )

simres = unnest( simres, cols = res )
simres
cli_alert_info( "Working with {nrow( simres )} rows of simulation data" )


table( simres$method )



##### Make wide-form for comparisons ######

names( simres )
reswide <- simres %>%
    dplyr::select( ID, runID, one_sided, N:scaled_C, method, LATE.hat ) %>%
    pivot_wider( names_from=method, values_from=LATE.hat )
reswide

dplyr::select( reswide, UNSTRAT:`2SLS` ) %>%
    cor( use = "pairwise.complete.obs" ) %>%
    round( digits=3 )


# How often to IV_a and IV_w agree?
mean( abs(reswide$IV_a - reswide$IV_w ) < 0.00000001 )
rm( reswide )



#### Calculate performance statistics ####


aggres <- simres %>%
    group_by( method, one_sided, scaled_C, het_tx) %>%
    summarise( mn = mean( LATE.hat, na.rm=TRUE),
               bias = mean( LATE.hat - LATE, na.rm=TRUE ),
               RMSE = sqrt( mean( (LATE.hat - LATE)^2, na.rm=TRUE ) ),
               SE = sd( LATE.hat, na.rm=TRUE ),
               n_c = mean( n_comp, na.rm=TRUE ),
               empty = 100*mean( empty > 0, na.rm=TRUE ), .groups = "drop" ) %>%
    group_by( scaled_C ) %>%
    mutate( per = SE / SE[method=="UNSTRAT"] ) %>%
    ungroup() %>%
    relocate( method, n_c )


skimr::skim( aggres )

table( is.na( aggres$n_c ), aggres$method )
table( is.na( aggres$empty ), aggres$method )
nrow(aggres)
filter( aggres, is.na( per ) )




#### Calculate R2 values #####

# Calculate predictive R2 of the covariate of compliance status  by
# generating giant datasets and running needed regressions and merge
# those into the performance data.

vals = simres %>%
    filter( het_tx == "yes" ) %>%
    dplyr::select( one_sided, pi_c:sd0 ) %>%
    unique()
vals
vals$N = 100000
nrow( vals )

# Check this aligns with the simulation code
pstrat = c( 0.35, 0.35, 0.25, 0.15 )
tau = c( 0, 0.2, 0.4, 0.6 )

dats <- pmap( vals, make_dat, pstrat = pstrat, tau=tau )
descs <- map_df( dats, poststratIV:::summarize_sim_data )
descs$scaled_C = vals$scaled_C
descs$one_sided = vals$one_sided

descs



# Note: We need to round so the merge works correctly.  Some minor
# precision issues were killing the merge before.
aggres$scaled_C = round( aggres$scaled_C, digits = 3 )
descs$scaled_C = round( descs$scaled_C, digits = 3 )
vals$scaled_C = round( vals$scaled_C, digits = 3 )

descs
aggres2 = left_join( aggres, descs, by=c( "scaled_C", "one_sided" ) )
aggres2



#### Look at expected bias of precision weighting #####


## Calculate the bias we would expect given precision weighting
length( dats )

biases <- map_df( dats, function( d ) {
    gg <- d %>%
        group_by( X ) %>%
        summarise( n = n() / 100,
                   p = mean(Z),
                   sigma2Y = var( Y0 ) / (n*p) + var( Y1 ) / (n*(1-p)),
                   f = mean( complier ),
                   wt = f^2 / sigma2Y,
                   LATE = mean( Y1 - Y0 ) / f )
    gg
    gg$LATEagg = mean( d$Y1 - d$Y0 ) / mean( d$S1 - d$S0 )
    gg$LATEwt = weighted.mean(gg$LATE, w = gg$wt )
    gg
}, .id = "group")

biases
table( biases$group )
length( dats )
filter( biases, X == "X2" )
filter( biases, group == 4 )
filter( biases, group == 1 )
filter( biases, group == 8 )


# Make table of biases
bias_table <- biases %>%
    dplyr::select( group, X, f, LATEagg, LATEwt ) %>%
    pivot_wider( names_from = "X", values_from= "f" ) %>%
    mutate( scaled_C = vals$scaled_C,
            one_sided = vals$one_sided,
            bias = round(LATEwt - LATEagg, digits=3 ) ) %>%
    rename( LATE = LATEagg ) %>%
    left_join( descs, by = c( "scaled_C", "one_sided" ) )
bias_table



map_df( dats, function( d ) {
    d %>% filter( complier == 1 ) %>%
        group_by( X ) %>%
        summarise( ETx = mean( Y1 - Y0 ),
                   sdTx = sd( Y1 - Y0 ),
                   sdY0 = sd( Y0 ),
                   n = n() ) %>%
        ungroup() %>%
        mutate( prop = n / sum(n),
                ATE = weighted.mean( ETx, n ) )
}, .id = "group" ) %>%
    knitr::kable()





##### Make plot of Bias, SE, RMSE #####


aggres2$one_sided = factor( aggres2$one_sided, c( TRUE, FALSE ), labels=c("one", "two" ) )
aggres2$method = fix_method_fct( aggres2$method )


aggresL <- aggres2 %>%
    pivot_longer( cols=c( "bias", "SE", "RMSE" ),
                  names_to = "metric",
                  values_to = "value" ) %>%
    mutate( metric = factor( metric, levels = c("bias", "SE", "RMSE" ) ) )

plt <- aggresL %>%
    filter( het_tx == "yes", method !="DSS" ) %>%
    ggplot( aes( R2_comp, value, col=method, lty=method, group=method ) ) +
    facet_grid( one_sided ~ metric ) +
    geom_hline( yintercept = 0 ) +
    geom_point() + geom_line() +
    labs( y = "Value",
          x = expression( paste( R^2, " of compliance") ),
          color = "", lty="" ) +
    theme_minimal() + theme( legend.position="bottom" )
add_color_line_thing(plt)


metric.labs = c( "Bias", "SE", "RMSE" )
names(metric.labs) = c("bias", "SE", "RMSE" )

aggresL = aggresL %>%
    mutate( method = fix_method_fct( method ) )

table( aggresL$method )


plt <- aggresL %>%
    filter( het_tx == "yes" ) %>%
    filter( method != "DSS" ) %>%
    ggplot( aes( R2_comp, value, col=method, lty=method ) ) +
    facet_grid( one_sided ~ metric, labeller = labeller(metric=metric.labs) ) +
    geom_hline( yintercept = 0 ) +
    geom_line() +
    labs( y = "Effect Size Units",
          x = expression( paste( R^2, " of compliance") ),
          color = "", lty="" ) +
    theme_minimal() +
    theme(plot.margin=margin(0,10,0,5),
          panel.spacing = unit(2, "lines") ) # Adjust space between facets


add_color_line_thing(plt)

ggsave( "figures/predC_Bias_SE_RMSE.pdf", width = 7, height = 3 )


