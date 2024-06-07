
# Analysis of the simulation focusing on violation of exclusion restriction.
#
# Simulation is run in simulation_exclusion_restriction.R


library( tidyverse )
library( scales )
library( poststratIV )
rm(list=ls())

theme_set( theme_minimal() )
source( "plot_helpers.R" )

library( cli )

cat( "\n\n" )
cli_alert_info( "Analyzing the simulation looking at violation of exclusion restriction\n" )


##### Load simulation data #####

simres = readRDS( here::here( "results/simulation_results_exclusion.rds" ) )

simres$err = NULL
simres$seed = NULL
simres$reps = NULL
simres$ID = 1:nrow(simres)

simres = unnest( simres, cols = res )
simres
nrow( simres )

table( simres$method )


#### Calculate performance statistics ####

if ( FALSE ) {
    skimr::skim(simres)
}

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


if ( FALSE ) {
    aggres %>%
        knitr::kable( digits=2) %>%
        print()
}



#### Calculate R2 values #####

# Calculate predictive R2 of the covariate of compliance status and
# merge those into the performance data

vals = simres %>%
    dplyr::select( one_sided, pi_c:sd0 ) %>%
    unique()
vals
vals$N = 100000
nrow( vals )

# Check this aligns with the simulation code

pstrat = c( 0.35, 0.35, 0.25, 0.15 )
tau = c( 0.5, 0.5, 0.5, 0.5 ) #0, 0.2, 0.4, 0.6 )

dats <- pmap( vals, make_dat, pstrat = pstrat, tau=tau )
descs <- map_df( dats, summarize_sim_data )
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


##### Make plot of Bias, SE, RMSE #####


aggres2 <- aggres2 %>%
    mutate( method = fct_reorder( method, per, mean, .desc=TRUE ) )# %>%
aggres2$one_sided = factor( aggres2$one_sided, c( TRUE, FALSE ), labels=c("one", "two" ) )


aggres2$method = fix_method_fct( aggres2$method )


aggresL <- aggres2 %>%
    pivot_longer( cols=c( "bias", "SE", "RMSE" ),
                  names_to = "metric",
                  values_to = "value" ) %>%
    mutate( metric = factor( metric, levels = c("bias", "SE", "RMSE" ) ) )


metric.labs = c( "Bias", "SE", "RMSE" )
names(metric.labs) = c("bias", "SE", "RMSE" )


plt <- aggresL %>%
    filter( method != "DSS" ) %>%
    ggplot( aes( R2_comp, value, col=method, lty=method ) ) + #, lty=het_tx) ) +
    # facet_wrap( ~ metric, nrow = 1 ) +
    facet_wrap( ~ metric, labeller = labeller(metric=metric.labs) ) +
    geom_hline( yintercept = 0 ) +
    #geom_point() +
    geom_line() +
    labs( y = "Effect Size Units",
          x = expression( paste( R^2, " of compliance") ),
          color = "Method", lty="Method" ) +
    theme_minimal() +
    #theme( legend.position="bottom" ) +
    #scale_x_continuous(labels = scales::percent) + # Format x-axis as percentages
    theme(plot.margin=margin(0,10,0,5),
          panel.spacing = unit(2, "lines") )# Adjust space between facets


plt <- add_color_line_thing(plt)
plt

  ggsave( "figures/exclusion_Bias_SE_RMSE.pdf", width = 7, height = 5 )


