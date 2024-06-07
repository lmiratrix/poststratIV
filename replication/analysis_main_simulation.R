
# This script processes the simulation results to generate several of
# the main plots of the paper.
#
# It merges both the two-sided and one-sided simulation scenarios into
# a single analysis.
#
# It can subset to either heterogeneous or not heterogeneous impact
# simulations, or just pool them together via pool.  See "SCENARIO"
# parameter at top.
#

library( tidyverse )
options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_minimal() )
library(RColorBrewer)

library( cli )

cat( "\n\n" )
cli_alert_info( "Analyzing the main simulation\n" )



rm(list = ls())

source( "plot_helpers.R" )

# Set whether to focus on subset of simulations, or just pool them
# all:
#SCENARIO = "nohet"
#SCENARIO = "het"
SCENARIO = "pool"



##### Utility functions #####


FIGURE_DIR = here::here( "figures/" ) #paste0( RESULT_DIR, "figures/" )
dir.create( FIGURE_DIR, showWarnings = FALSE )
#  FIGURE_DIR = here::here( "../JRSS-A submission/unblinded - JRSS-A with reviewer edits/figures/" )

cat( "Saving figures to ", FIGURE_DIR, "\n" )

make_file = function( name,
                      path = FIGURE_DIR ) {
    paste0( path, name, "_", SCENARIO, "_dual.pdf" )
}




#### Load Simulation results  ####

simres = readRDS( "results_onesided/simulation_results.rds" )
simres2 = readRDS( "results_twosided/simulation_results.rds" )
simres = bind_rows( one = simres, two = simres2, .id = "sides" )
rm( simres2 )



# Examine errors
errs = map_dbl( simres$err, length )
simres$err[[1]]
table( errs )

simres$err = NULL
simres$seed = NULL
simres$reps = NULL

simres %>%
    dplyr::select( -chunkNo, -res ) %>%
    unique()


simres = unnest( simres, cols = res )
simres
cat( "# rows of results:", nrow( simres ), "\n" )
og_nrow = nrow( simres )


# Patch very small pi.hat giving crazy LATE.hat estimates
if ( TRUE ) {
    #simres$pi.hat = zapsmall( simres$pi.hat, digits = 4 )
    pp = abs( simres$pi.hat ) < 0.0001
    table(pp)
    ss = !is.na( simres$LATE.hat )
    table(ss)
    simres$LATE.hat[ pp ] = NA
    ss = !is.na( simres$LATE.hat )
    table(ss)
}


miss <- with( simres, table( method, is.na( LATE.hat ) ) )
t = apply( miss, 1, sum )
round( 100 * miss / t, digits = 2 )

simres = filter( simres, !is.na( LATE.hat ) )
simres



simres <- dplyr::select( simres, -SE.ITT, -cov.ITT.pi, -SE.pi,
                         -ITT.hat.ps, -pi.hat.ps, -SE.ITT.ps, -SE.pi.ps,
                         -cov.ITT.pi.ps )
simres

cat( "# rows of results (final):", nrow( simres ), "\n",
     "\t", 100*round( 1 - nrow(simres) / og_nrow, digits = 2 ), "percent dropped\n" )


#### Handle extreme value estimates ####

if ( FALSE ) {
    # Look at extreme values--how extreme are they?
    filter( simres, LATE.hat > 40 ) %>%
        dplyr::select( sides:het_tx, one_sided, method, pi.hat, LATE.hat ) %>%
        arrange( -LATE.hat )
}


if ( TRUE ) {

    # Max out LATE estimates at THRESHOLD
    THRESHOLD = 10
    simres = mutate( simres,
                     LATE.hat = pmax( -THRESHOLD, pmin( LATE.hat, THRESHOLD ) ),
                     SE.wald = pmin( abs( SE.wald ), THRESHOLD ),
                     SE.delta = pmin( abs( SE.delta ), THRESHOLD ),
                     maxed = abs(LATE.hat) == THRESHOLD,
                     maxed_SE = (abs(SE.wald) == THRESHOLD) | (abs(SE.delta) == THRESHOLD) )

    table( maxed=simres$maxed, isna = is.na(simres$LATE.hat), useNA="always" )
    table( PEmax=simres$maxed, SEmax=simres$maxed_SE, useNA="always" )
    mean( simres$maxed_SE, na.rm=TRUE )
    table( is.na( simres$SE.delta ), simres$method )
    table( is.na( simres$SE.wald ), simres$method )
    simres$SE.delta[ simres$maxed ] = NA
    simres$SE.wald[ simres$maxed ] = NA

    # length( unique( simres$SE.delta ) )
    # length( simres$SE.delta )
    # tt = table( simres$SE.delta )
    # tt[ tt > 5 ]
    # head( tt[ tt > 4 ] )
    # mean( tt == 1 )
    # table( table( simres$SE.delta ) )

    #simres$SE.delta
    #head( sort( simres$SE.delta, decreasing = TRUE ) )
    #mean( simres$SE.delta == 40, na.rm = TRUE )
    #mean( simres$SE.delta > 40, na.rm = TRUE )

    summary( simres$SE.delta )

    # overall percent of extreme estimates
    mean( simres$maxed )
    mean( simres$maxed_SE, na.rm=TRUE )

    # Percent of values at threshold
    simres %>% group_by( method, N, sides ) %>%
        summarise( maxed = 100 * mean( maxed ) ) %>%
        pivot_wider( names_from=c( "N", "sides" ), values_from="maxed" ) %>%
        knitr::kable( digits = 1 )

    simres %>% group_by( method, sides ) %>%
        summarise( maxed = 100 * mean( maxed ),
                   maxSE = 100 * mean( maxed_SE )) %>%
        pivot_wider( names_from=c( "sides" ), values_from=c( "maxed", "maxSE" ) ) %>%
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


    mm <- simres %>%
        group_by( method, N, pi_c ) %>%
        summarise( max_est = 100 * mean( maxed, na.rm=TRUE ),
                   max_SE = 100 * mean( maxed_SE, na.rm = TRUE ),
                   max = 100 * mean( maxed | maxed_SE, na.rm=TRUE )   )
    #pivot_wider( names_from="pi_c", values_from=c( "max", "max_SE", "max_est" ) )

    table( maxPE = simres$maxed, maxSE = simres$maxed_SE )

    mm %>%
        pivot_wider( names_from="pi_c", values_from=c(  "max_SE", "max" ) ) %>%
        filter( method != "Oracle" ) %>%
        arrange( N, method ) %>%
        knitr::kable( digits= 1 )
    mm
    mm %>%
        dplyr::filter( method != "Oracle" ) %>%
        pivot_longer( cols = c( max, max_est, max_SE ) ) %>%
        ggplot( aes( pi_c, value, col = as.factor(N) ) ) +
        facet_grid( method ~ name ) +
        geom_line()


    simres %>% group_by( sides, N, pi_c, nt_shift, pred_comp, one_sided, method ) %>%
        summarise( maxPE = 100 * mean( maxed, na.rm=TRUE ),
                   maxSE = 100 * mean( maxed_SE, na.rm=TRUE ),
                   .groups = "drop" ) %>%
        group_by( method ) %>%
        summarise( meanPE = mean(maxPE),
                   maxPE = max(maxPE),
                   meanSE = mean(maxSE),
                   maxSE = max(maxSE) )
}


if ( FALSE ) {
    # Peek at distribution of impact estimates

    simres %>%
        filter( method == "UNSTRAT" ) %>%
        ggplot( aes( LATE.hat - LATE ) ) +
        facet_grid( pred_comp + pred_Y ~ N ) +
        geom_histogram() +
        geom_vline( xintercept = 0, col="red")

}

table( simres$method )

# counts of each simulation scenario type
counts <- simres %>%
    group_by( N, pi_c, nt_shift, pred_comp, pred_Y, het_tx, sd0, sides, method ) %>%
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



#skimr::skim(simres)
filter( simres, is.na( SE.delta ) ) %>%
    dplyr::select( -(N:sd0) )

mean( is.na( simres$SE.wald ) )


summary( simres$SE.delta - simres$SE.wald )


res <- simres %>%
    group_by( method, N, pi_c, nt_shift, sd0, pred_comp, pred_Y, het_tx, sides ) %>%
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


summary( res$SE^2 / res$E_SE2_delta)



##### Very high level stats across everything #####


aggstat <- res %>%
    group_by( N, pi_c, method, sides ) %>%
    summarise( SE = mean( SE, na.rm=TRUE ),
               abs_bias = mean( abs_bias ),
               RMSE = mean( RMSE ),
               SEwald = sqrt( mean( E_SE2_wald ) ),
               SEdelta = sqrt( mean( E_SE2_delta ) ),
               sdwald = mean( sd_SE_wald ),
               sddelta = mean( sd_SE_delta ),
               .groups = "drop" )

aggstat = mutate( aggstat, Nc = N * pi_c )

ss <- aggstat %>%
    dplyr::select( N, pi_c, Nc, method, SE, sides ) %>%
    pivot_wider( names_from="method", values_from="SE" ) %>%
    arrange( sides, Nc ) %>%
    mutate( ratio = UNSTRAT/Oracle )
ss


if ( FALSE ) {
    # How many times larger is the unstratified's true SE vs. ideal of the
    # oracle (perfect identification of compliers)
    ggplot( ss, aes( pi_c, ratio, col=as.factor(N) ) ) +
        geom_point( size = 3) +
        geom_hline(yintercept = 1 )



    ggplot( aggstat, aes( Nc, SE, col = as.factor( N ) ) ) +
        facet_wrap( ~ method ) +
        geom_hline( yintercept = 1, col="grey" ) +
        geom_jitter( width = 3.5, size=3 ) +
        coord_cartesian( ylim=c(0,2.5 ) )
}



##### Look at triple of bias, se, RMSE of estimators ######

rlong = res %>%
    dplyr::select( method, N, pi_c, nt_shift, sd0, abs_bias,
                   SE, RMSE, pred_comp, pred_Y, het_tx, sides ) %>%
    pivot_longer( cols=c(abs_bias, SE, RMSE),
                  names_to="metric", values_to="value" )
rlong
rlong = mutate( rlong,
                metric = factor( metric, levels=c("abs_bias", "SE", "RMSE" ) ) )




# What scenarios totally fail?
filter( rlong, value > 200 )


# summary of overall performance
rlong
stats <- rlong %>% group_by( method, metric, sides, N ) %>%
    summarise( Emetric = mean(value),
               sdmetric = sd(value) )
stats


metric.labs = c( "|Bias|", "SE", "RMSE" )
names(metric.labs) = c("abs_bias", "SE", "RMSE" )



##### Overall results figure in paper #####
rlong = mutate( rlong,
                category = fct_collapse( method,
                                         baseline = c( "UNSTRAT", "Oracle" ),
                                         core = c( "IV_a", "IV_w", "2SLS" ),
                                         weight = c( "DSS0", "PWIV", "DSF") )
)

rlong$method = fix_method_fct( rlong$method )

plt <- rlong %>%
    mutate(
        #        method = fct_reorder2( method, N * (metric=="RMSE"), value, function( x, y ) {
        #      a = y[ x == max(x) ]
        #      mean(a) }, .desc=TRUE ),
        N = as.factor( N ) ) %>%
    ggplot( aes( N, value, col=method, lty=method, group=method ) ) +
    facet_grid( sides ~ metric, labeller = labeller(metric=metric.labs) ) +
    #facet_grid( . ~ metric + sides, labeller = labeller(metric=metric.labs) ) +
    labs( y = "Effect Size Units", color = "Method", lty="Method" ) +
    scale_x_discrete( breaks = c( 500, 1000, 2000 ), labels = c( "500", "1K", "2K" ),
                      expand = c( 0, 0.1 ) ) +
    stat_summary(fun=mean, geom="line") +
    theme(panel.spacing = unit(1, "lines"),
          panel.border = element_rect(colour = "grey", fill=NA, linewidth=1) )

plt <- add_color_line_thing( plt )
plt

ggsave( file=make_file("overall_performance_plot"), width = 7, height = 3 )



#### Bias analysis ####


rlong_agg <- rlong %>%
    filter( het_tx == "yes" ) %>%
    group_by( method, metric, N, pi_c, sides ) %>%
    summarize( value = mean( value ) ) %>%
    group_by( metric, N, pi_c ) %>%
    mutate( perval = value / value[method=="Unstrat"] )

rlong_agg %>%
    mutate( method = fix_method_fct( method ),
            method = fct_relevel( method, "Unstrat", "IV[a]", "`2SLS`", "IV[w]", "DSS0" )) %>%
    filter( metric == "abs_bias",
            method != "Oracle" ) %>% # method != "Oracle", method != "UNSTRAT" ) %>%
    ggplot( aes( pi_c, value,  col=as.factor(N) ) ) + #,  col=as.factor(nt_shift) ) ) +
    facet_grid( sides ~ method, labeller = label_parsed ) +
    geom_hline( yintercept = 0 ) +
    geom_line() + geom_point() +
    # coord_cartesian( xlim = c( 0.04, 0.11 ), ylim=c(0,2.5) ) +
    scale_x_continuous( breaks = unique( rlong_agg$pi_c ), labels = scales::percent ) +
    labs( x = expression( pi[c] ), y = "Bias", col = "N" ) +
    #scale_y_continuous( labels = scales::percent ) +
    theme(panel.spacing = unit(1, "lines")) +
    #  geom_hline( yintercept = 1, lty = 2) +
    theme( legend.position="bottom" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_brewer(palette = "Set2") # Use 'Set2' palette

ggsave( file=make_file( "bias_ntshift_plot" ), width = 7, height = 3 )


#### Looking at SE-hat performance: calibration #####

cat( "Now looking at SE-hat's performance\n" )


res

summary( res$SE )

# Analysis: Are standard errors too big or too small?
calib_calc <- res %>%
    filter( method != "Oracle" ) %>%
    group_by( method, N ) %>%
    dplyr::select( method, N, E_SE2_delta, E_SE2_wald, med_SE_delta, med_SE_wald, SE, sides) %>%
    pivot_longer( col=c( E_SE2_delta, E_SE2_wald, med_SE_delta, med_SE_wald ),
                  names_to = c(".value", "estimator"), names_pattern = '(.*_.*)_(.*)' ) %>%
    #names_prefix="E_SE2_", names_to="estimator", values_to="E_SE2" ) %>%
    mutate( SE2_ratio = E_SE2 / SE^2,
            med_ratio = med_SE / SE,
            SE_ratio = sqrt( SE2_ratio ) ) %>%
    ungroup()


calib_calc

# Drop double-entries of 2SLS estimators
calib_calc = calib_calc %>%
    filter( estimator == "wald" | method != "2SLS" ) %>%
    mutate( estimator = ifelse( method == "2SLS", "2SLS", estimator ) )
table( calib_calc$method, calib_calc$estimator )


calib_calc$method = fix_method_fct( calib_calc$method )
calib_calc %>%
    mutate( estimator = ifelse( estimator=="wald", "Bloom", "Delta" ) ) %>%
    mutate( estimator = ifelse( method == "2SLS", "2SLS", estimator ) ) %>%
    ggplot( aes( as.factor(N), SE_ratio, col=estimator ) ) +
    facet_grid( sides ~ method , scales="free", labeller = label_parsed ) +
    geom_hline( yintercept = 1, col="grey" )  +
    geom_boxplot( width = 0.5) +
    labs( x = "Sample Size (N)", y = "Calibration (%)", col="Estimator" ) +
    scale_y_continuous( labels = scales::percent_format( ) ) +
    scale_x_discrete( breaks = c( 500, 1000, 2000 ), labels = c( "500", "1K", "2K" ) ) +
    scale_color_brewer(palette = "Set2") # Use 'Set2' palette

ggsave( file=make_file("se_estimator_plot" ), width = 7, height = 3 )


# Comparing delta vs wald
calib_calc
aa <- calib_calc %>%
    dplyr::select( method, N, sides, estimator, SE, E_SE2 ) %>%
    pivot_wider(names_from=estimator, values_from=E_SE2 )



nrow(aa)
table( aa$method )
sample_n( aa, 20 ) %>% arrange( method )
summary( aa$delta - aa$wald )



##### Looking at median estimated SE #####

calib_calc

#calib_calc = mutate( calib_calc,
#                     med_ratio = med_ratio - 1 )


calib_calc$method = fix_method_fct( calib_calc$method )

calib_calc %>%
    mutate( estimator = ifelse( estimator=="wald", "Bloom", "Delta" ) ) %>%
    mutate( estimator = ifelse( method == "2SLS", "2SLS", estimator ) ) %>%
    ggplot( aes( as.factor(N), med_ratio, col=estimator ) ) +
    facet_grid( sides ~ method , scales="free", labeller = label_parsed ) +
    geom_hline( yintercept = 1, col="grey" )  +
    geom_boxplot( width = 0.5) +
    labs( x = "Sample Size (N)", y = "Calibration (%)", col="Estimator" ) +
    scale_y_continuous( labels = scales::percent_format( ) ) +
    scale_x_discrete( breaks = c( 500, 1000, 2000 ), labels = c( "500", "1K", "2K" ) ) +
    scale_color_brewer(palette = "Set2") # Use 'Set2' palette

ggsave( file=make_file("se_estimator_plot_med" ), width = 7, height = 3 )


arrange( calib_calc, med_ratio )

# Distribution of SE estimates across simulations is skewed.
gg = filter( simres, N == 500 )
gg = gg %>%
    pivot_longer( cols=c(SE.wald, SE.delta),
                  names_to = "SE_method", values_to="SE.hat" ) %>%
    mutate( SE.hat2 = SE.hat^2 )

ggplot( gg, aes( SE.hat2 ) ) +
    facet_grid( SE_method ~ pred_Y + pred_comp + het_tx, labeller = label_both )  +
    geom_histogram() +
    scale_x_continuous( limits=c(0,10 ) )

gg %>% group_by( pred_Y, pred_comp, het_tx, SE_method ) %>%
    summarise( sdSE2 = sd( SE.hat2, na.rm=TRUE ) ) %>%
    pivot_wider( names_from=SE_method, values_from=sdSE2 )


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
    facet_grid( sides ~ SEmethod ) +
    geom_boxplot()




a = filter( resSE, N == 500, pi_c == 0.05, nt_shift == 0, sd0 == 1, het_tx =="no",
            pred_comp=="no", pred_Y =="yes", SEmethod=="wald", sides =="two" )
a
stopifnot( sum( duplicated( a$method ) ) == 0 )


ratios <- resSE %>%
    group_by(  N, pi_c, nt_shift, sd0, pred_comp, pred_Y, het_tx, SEmethod, sides ) %>%
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






##### Comparing standard errors to unstratified baseline #####

cat( "Comparing SE as ratio to baseline\n" )

ratios <- simres %>%
    group_by(  method, N, pi_c, nt_shift, pred_comp, pred_Y, het_tx, sd0, sides ) %>%
    summarise( var_f = var( pi.hat ),
               var_ITT = var( ITT.hat ),
               var_LATE = var( LATE.hat ), .groups = "drop" ) %>%
    pivot_longer( cols=var_f:var_LATE,
                  names_to = "estimator",
                  names_prefix="var_",
                  values_to="var" ) %>%
    group_by( N, pi_c, nt_shift, pred_comp, pred_Y, sd0, estimator ) %>%
    mutate( ratio = var / var[method=="UNSTRAT"],
            n = n() ) %>%
    ungroup()
ratios

filter( ratios, N==500, nt_shift==-1, pred_comp=="yes", pi_c == 0.05,
        pred_Y == "no" ) %>%
    arrange( estimator )

ratios = filter( ratios, method != "UNSTRAT" )
ggplot( sample_n( ratios, nrow(ratios) ),
        aes( pred_comp, ratio, group = pred_Y, col=pred_Y ) ) +
    facet_grid( method ~ estimator ) +
    geom_jitter( width=0.1 )


var_rat = filter( ratios, estimator == "LATE" )
ggplot( var_rat, aes( pi_c, ratio, col = pred_Y, pch=pred_comp ) ) +
    facet_wrap( ~ method ) +
    geom_jitter( width=0.001 )


# Which scenarios lead to superior estimation of the compliance rate?
gg = filter( ratios, estimator == "f", method=="IV_a" )
gg
gg %>% group_by( pred_comp, pred_Y, pi_c ) %>%
    summarize( ratio = mean( ratio ) ) %>%
    arrange( pi_c )
gg %>% group_by( pred_comp, pred_Y ) %>%
    summarize( ratio = mean( ratio ) )


ratio_w = filter( ratios, estimator == "LATE" )

filter( ratio_w, N == 500, pi_c == 0.05, pred_comp == "yes" )

ratio_w_agg <- ratio_w %>%
    group_by( method, N, pi_c, pred_comp, pred_Y ) %>%
    summarise( ratio = mean(ratio) )


library( scales )

ratio_w_agg %>%
    filter( method %in% c( "IV_w", "IV_a" ) ) %>%
    mutate( pred_comp = str_to_title(pred_comp),
            pred_Y = str_to_title(pred_Y),
            #  pred_Y = ifelse( pred_Y == "yes",
            #                   "Y predictive", "Y not predictive" ),
            method = recode_factor(method,
                                   `IV_w` = "IV[w]",
                                   `IV_a` = "IV[a]" ) ) %>%
    ggplot( aes( pi_c, ratio, col = pred_comp, lty = pred_Y ) ) +
    facet_grid( method ~ N, labeller = label_parsed ) + # label_bquote(rows = .(method)) ) +
    geom_point() + geom_line() +
    #theme( panel.border = element_rect(color = "black", fill = NA, linewidth = 1) ) +
    coord_cartesian( ylim=c(0,1) ) +
    labs( x = "Proportion complier",
          #y = expression( paste( "% Variance of ", IV[w], " vs. Unstratified" ) ),
          y = "% Variance vs. Unstratified",
          lty = "Predict Y",
          col = "Predict C" ) +
    scale_x_continuous(  breaks=unique( ratio_w$pi_c ),
                         labels = scales::percent_format( ) ) +
    scale_y_continuous( labels = scales::percent_format( accuracy = 1 ) ) +
    #scale_linetype(labels = scales::parse_format()) +
    theme(panel.spacing = unit(2, "lines"))

ggsave( file=make_file( "variance_ratio_plot_v2" ), width = 8, height = 3 )




