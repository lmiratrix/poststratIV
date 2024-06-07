
# This script makes the SE instability plot in the appendix of the paper
#
# It can analyze either het or nohet simulations, or just pool them
# together via pool.  See "SCENARIO" parameter at top.
#


library( tidyverse )
library( cli )
options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_minimal() )


# Set TRUE for one-sided simulation
ONE_SIDED_SIMULATION = TRUE



#SCENARIO = "nohet"
#SCENARIO = "het"
SCENARIO = "pool"


cat( "\n\n" )
cli_alert_info( "Analyzing SE instability with one-sided={ONE_SIDED_SIMULATION} and SCENARIO={SCENARIO}\n" )



RESULT_DIR = here::here( "results_twosided/" )
if ( ONE_SIDED_SIMULATION ) {
    RESULT_DIR = here::here( "results_onesided/" )
}


##### Utility functions #####


FIGURE_DIR = here::here( "figures/" ) #paste0( RESULT_DIR, "figures/" )
dir.create( FIGURE_DIR, showWarnings = FALSE )
#  FIGURE_DIR = here::here( "../JRSS-A submission/unblinded - JRSS-A with reviewer edits/figures/" )

cat( "Saving figures to ", FIGURE_DIR, "\n" )

make_file = function( name,
                      path = FIGURE_DIR ) {
    add = "twosided"
    if ( ONE_SIDED_SIMULATION ) {
        add = "onesided"
    }
    paste0( path, name, "_", SCENARIO, "_", add, ".pdf" )
}




#### Load Simulation results  ####

# Results loaded from
simres = readRDS( paste0( RESULT_DIR, "simulation_results.rds" ) )
simres

# Examine errors
errs = map_dbl( simres$err, length )
simres$err[[1]]
table( errs )

simres$err = NULL
simres$seed = NULL
simres$reps = NULL

simres = unnest( simres, cols = res )
simres
cat( "# rows of results:", nrow( simres ), "\n" )


with( simres, table( method, is.na( LATE.hat ) ) )
simres = filter( simres, !is.na( LATE.hat ) )
simres



simres <- dplyr::select( simres, -SE.ITT, -cov.ITT.pi, -SE.pi,
                         -ITT.hat.ps, -pi.hat.ps, -SE.ITT.ps, -SE.pi.ps,
                         -cov.ITT.pi.ps )
simres

cat( "# rows of results (final):", nrow( simres ), "\n" )


#### Handle extreme value estimates ####
if ( TRUE ) {

    # Max out LATE estimates at THRESHOLD
    THRESHOLD = 10
    simres = mutate( simres,
                     LATE.hat = pmax( -THRESHOLD, pmin( LATE.hat, THRESHOLD ) ),
                     maxed = abs(LATE.hat) == THRESHOLD )


    # overall percent of extreme estimates
    mean( simres$maxed )

    # Percent of values at threshold
    simres %>% group_by( method, N ) %>%
        summarise( maxed = 100 * mean( maxed ) ) %>%
        pivot_wider( names_from="N", values_from="maxed" ) %>%
        knitr::kable( digits = 1 )



    # Look at distribution of SE estimates
    simres %>% group_by( method, pi_c, N ) %>%
        summarise( propneg = mean( pi.hat <= 0, na.rm=TRUE ),
                   min_SE = min( SE.wald, na.rm=TRUE ),
                   max_SE = max( SE.wald, na.rm=TRUE )) %>%
        arrange( -propneg )

    filter( simres, pi.hat < 0 ) %>%
        dplyr::select( N:method, pi.hat ) %>%
        sample_n( n() )

    simres <- simres %>% mutate(
        SE.wald = pmax( -2*THRESHOLD, pmin( abs( SE.wald ), 4*THRESHOLD ) ),
        SE.delta = pmax( -2*THRESHOLD, pmin( abs( SE.delta ), 4*THRESHOLD ) ),
        maxed_SE = abs(SE.wald) == 2*THRESHOLD | abs(SE.delta) == 4 * THRESHOLD )

    simres %>% group_by( method, N, pi_c ) %>%
        summarise( max_est = 100 * mean( maxed, na.rm=TRUE ),
                   max_SE = 100 * mean( maxed_SE, na.rm = TRUE ),
                   max = 100 * mean( maxed | maxed_SE, na.rm=TRUE )   )%>%
        pivot_wider( names_from="pi_c", values_from=c( "max", "max_SE", "max_est" ) ) %>%
        filter( method != "Oracle" ) %>%
        arrange( N, method ) %>%
        knitr::kable( digits= 1 )
}



table( simres$method )

# counts of each simulation scenario type
counts <- simres %>%
    group_by( N, pi_c, nt_shift, pred_comp, pred_Y, het_tx, sd0, method ) %>%
    summarise( n = n() )

counts %>%
    group_by( method ) %>%
    summarise( n_min = min( n ),
               n_max = max( n ),
               nbar = mean( n ) )
filter( counts, n <= 500 )





##### Subset to desired scenario collection ######

#simres = filter( simres, method != "Oracle" )
if ( SCENARIO == "nohet" ) {
    cat( "Subsetting to het_tx == 'no'\n" )
    simres = filter( simres, het_tx == "no" )
} else if ( SCENARIO == "het" ) {
    cat( "Subsetting to het_tx == 'yes'\n" )
    simres = filter( simres, het_tx == "yes" )
} else {
    cat( "Not subsetting, will end up pooling het and no het sims\n" )
}

table( simres$method )
simres = simres %>%
    filter( method != "DSS" ) %>%
    mutate( method = factor( method,
                             levels = c( "UNSTRAT",
                                         "IV_w", "IV_a",
                                         "DSS0",
                                         "2SLS",
                                         "DSS", "DSF", "PWIV",
                                         "Oracle" ) ) )






#### Make aggregate performance metrics  #####

filter( simres, is.na( SE.delta ) ) %>%
    dplyr::select( -(N:sd0) )

mean( is.na( simres$SE.wald ) )


summary( simres$SE.delta - simres$SE.wald )


res <- simres %>% group_by( method, N, pi_c, nt_shift, sd0, pred_comp, pred_Y, het_tx ) %>%
    filter( !is.na( SE.wald ) ) %>%
    summarise( E_est = mean( LATE.hat ),
               bias = mean( E_est - LATE ),
               abs_bias = abs( bias ),
               SE = sd( LATE.hat ),
               RMSE = sqrt( mean( (LATE.hat-LATE)^2 ) ),
               E_SE_wald = mean( SE.wald ),
               E_SE2_wald = mean( SE.wald^2 ),
               med_SE_wald = median( SE.wald ),
               E_SE_delta = mean( SE.delta ),
               E_SE2_delta = mean( SE.delta^2 ),
               med_SE_delta = median( SE.delta ),
               sd_SE_wald = sd( SE.wald ),
               sd_SE_delta = sd( SE.delta ),
               .groups = "drop" )

res

res = mutate( res, method = factor( method ) ) %>%
    mutate( method = relevel( method, ref = "UNSTRAT" ) )





##### Looking at variability of SEhat estimator ######

resSE <- res %>%
    filter( method != "Oracle" ) %>%
    pivot_longer( cols=E_SE_wald:sd_SE_delta,
                  names_pattern = c("(.*)_(.*)"),
                  names_to=c( ".value", "SEmethod" ) ) %>%
    dplyr::select( -E_est, -bias, -abs_bias, -RMSE )

resSE

# Standard deviation of SE-hat as percent of true SE
resSE <- resSE %>%
    mutate( stability = 100*sd_SE / SE )
resSE

ggplot( resSE, aes( method, stability ) ) +
    facet_wrap( ~ SEmethod ) +
    geom_boxplot()





a = filter( resSE, N == 500, pi_c == 0.05, nt_shift == 0, sd0 == 1, het_tx =="no",
            pred_comp=="no", pred_Y =="yes", SEmethod=="wald")
a
stopifnot( sum( duplicated( a$method ) ) == 0 )


ratios <- resSE %>%
    group_by(  N, pi_c, nt_shift, sd0, pred_comp, pred_Y, het_tx, SEmethod ) %>%
    mutate( dir_ratio = sd_SE / sd_SE[method=="UNSTRAT"],
            ratio = stability / stability[method=="UNSTRAT"],
            n = n() ) %>%
    ungroup()
sample_n( ratios, 10 )

filter( ratios, N==500, nt_shift==-1, pred_comp=="yes", pi_c == 0.05,
        pred_Y == "no" ) %>%
    arrange( method )

ratios = filter( ratios, method != "UNSTRAT" ) %>%
    mutate( ratio = ratio - 1,
            SEmethod = ifelse( SEmethod=="wald", "Bloom", "Delta" ) )


ggplot( ratios,
        aes( as.factor(N), ratio, col=SEmethod  ) ) +
    #facet_grid( method ~ SEmethod ) +
    facet_wrap( ~ method, nrow=1 ) +
    #geom_jitter( width=0.1 ) +
    geom_boxplot( width = 0.5 ) +
    geom_hline( yintercept = 0 ) +
    labs( x = "Sample Size", y = "Relative stability",
          col = "Method" ) + #y = "[sd_SE_p / SE_p] / [sd_SE_u / SE_u]") +
    scale_y_continuous( labels = scales::percent_format( ) ) +
    scale_x_discrete( breaks = c( 500, 1000, 2000 ), labels = c( "500", "1K", "2K" ),
                      expand = c( 0, 0.1 ) ) +
    theme(plot.margin=margin(0,10,0,5),
          panel.spacing = unit(2, "lines") ) # Adjust space between facets



ggsave( file=make_file("se_instability_plot" ), width = 7, height = 3 )







