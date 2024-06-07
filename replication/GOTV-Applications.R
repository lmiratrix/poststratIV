
# Code for data illustrations in "Improving instrumental variable estimators with post-stratification"
# 2024, by Pashley, Keele, & Miratrix (unless otherwise noted)

library( foreign )
library( xtable )
library( tidyverse )
library ( AER )
library( here )
library( Matrix )
library( lme4 )

# Load package for post-stratified analysis
library( poststratIV )

rm(list=ls())



#### Analyze the GOTV experiment of Gerber & Green 2003 (door-to-door canvassing) ####

# Data can be found at: https://isps.yale.edu/research/data/d017
# Put data in local folder GOTV


##### Load and prep data #####
data <- read.dta( here::here( "GOTV/local.dta") )
data <- subset(data, wellbehaved==1)

##### Conduct analysis #####
# Stratify on the "turf" variable
# Summarize data based on "turf"
table.desc <- data %>% group_by( turf ) %>%
    mutate( n = n() ) %>%
    group_by( turf, treatmen ) %>%
  summarise( n = n[[1]],
             p_con = mean( contact ),
             p_vote = mean( voted01 ), .groups="drop" ) %>%
  pivot_wider( names_from=treatmen, values_from=c( p_con, p_vote ) ) %>%
    mutate( complier = p_con_1 - p_con_0 )

table.desc

summary(table.desc$n)
summary( table.desc$p_vote_0 )

# For descriptive purposes, examine the range of voting prevalence across turfs
dd <- data %>% group_by( turf )
summary( data$voted01 )
M = glmer( voted01 ~ 1 + (1|turf),
           family = binomial,
           data = filter( data, treatmen == 0 ) )
arm::display( M )
re = ranef( M )$turf
dd = data %>% group_by( turf ) %>%
    summarise( n = n(),
               pvote = mean( voted01 ),
               ptx = mean( treatmen ) )

dd$pds = predict( M, newdata=dd, type="response" )
ggplot( dd, aes( pvote, pds, size=n ) ) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1 )
summary( dd$pds )
table(data$treatmen)


# Estimate impacts
iv.s <- IV.est.strat(data, Yobs = "voted01", Z="treatmen", S="contact",
                        strat_var = "turf",
                     tidy_table = TRUE )

# Reformat results to include test statistics, p-value, and SE reduction
iv.s <- iv.s %>%
    mutate(t = LATE.hat / SE.wald,
           pvalue = 2 * pmin( 1 - pnorm(t), pnorm( t )),
           rel_SE = 100 * SE.wald / SE.wald[ Xblk=="UNSTRAT"] )

# Filter to LATE results (i.e., exclude individual turf results)
iv_key <- iv.s  %>%
    dplyr::select( Xblk, pi.hat, LATE.hat, SE.wald, rel_SE, n, pvalue ) %>%
    filter(row.names(iv.s) %in% c("151", "152", "153", "154", "155", "156", "157"))


iv_key

# Format into table
print(xtable(iv_key, digits = c(0,0,2,3,3,0, 0,2 )), include.rownames=FALSE )

# Look at summary statistics for turfs
summary(iv.s[1:150,3])




#### Analyze the GOTV experiment of Gerber & Green 2000 (door-to-door canvassing) ####

##### Load and prep data #####

# Data and code for prepping from Hansen & Bowers (2009)
# available at
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/Q6CID7
hh_dat <- read.dta(here::here( "GOTV/NHrep_household.dta"),convert.underscore=TRUE)
ivdl_dat <- read.dta( here::here( "GOTV/NHrep_individual.dta" ),convert.underscore=TRUE)

hh_dat$AGE1[hh_dat$AGE1MISS==1] <- NA
hh_dat$AGE2[hh_dat$AGE2MISS==1] <- NA
hh_dat$V98.1[hh_dat$V98.1==99] <- NA
hh_dat$V98.2[hh_dat$V98.2==99] <- NA
hh_dat$V96.2.0[hh_dat$PERSONS==1] <- hh_dat$V96.1.0[hh_dat$PERSONS==1]
hh_dat$hh.n <-   hh_dat$PERSONS
hh_dat$hh_turnout <- ifelse(hh_dat$PERSONS==1, ifelse(is.na(hh_dat$V98.1),0,hh_dat$V98.1),
                            ifelse(is.na(hh_dat$V98.1),0,hh_dat$V98.1) +
                                ifelse(is.na(hh_dat$V98.2),0,hh_dat$V98.2) )
hh_dat$ipc <- as.logical(ivdl_dat$cntany[match(hh_dat$ID1, ivdl_dat$id1)])
hh_dat$phc <- as.logical(ivdl_dat$pcntany[match(hh_dat$ID1, ivdl_dat$id1)])
hh_dat$mailings <- ivdl_dat$mailings[match(hh_dat$ID1, ivdl_dat$id1)]
hh_dat$phongotv <- as.logical(hh_dat$PHONGOTV)
hh_dat$persngrp <- as.logical(hh_dat$PERSNGRP)
hh_dat$mailgrp <- as.logical(hh_dat$MAILGRP)
row.names(hh_dat) <- as.character(hh_dat$ID1)
hh_dat$ccs <- complete.cases(hh_dat$hh_turnout)
hh_dat$phonegrp <- as.numeric(hh_dat$PHONGOTV | hh_dat$BLOOD!="not selected")
hh_dat$any_to <- as.numeric(hh_dat$hh_turnout > 0)
hh_dat$any_to_96 <- as.numeric(hh_dat$V96.1.1 + hh_dat$V96.2.1 > 0)


##### Conduct analysis #####

# Stratify based on prior vote behavior, household size, and age
hh_dat2 <- hh_dat
# Initialize stratification variable based on
# prior vote behavior and household size
hh_dat2$blk <- interaction(hh_dat2$any_to_96, hh_dat2$PERSONS)

# Discretize (average household) age to use in stratification
hh_dat2 <- hh_dat2  %>%  mutate(agerank = ntile(0.5 * (AGE1 + AGE2),4))
# Updated stratification variable to include age
hh_dat2$blk2 <-interaction(hh_dat2$blk, hh_dat2$agerank)

# Remove observations with relevant missing data for simplicity
hh_dat2 <- hh_dat2[complete.cases(hh_dat2[,c("any_to", "ipc", "persngrp")]),]

# Estimate impacts (ignoring other factors)
strata_blk2_results1 <- IV.est.strat( data = hh_dat2,
                                      Yobs="any_to",
                                      S = "ipc", Z="persngrp", strat_var = "blk2" )
dim( strata_blk2_results1 )

# Print out core results
strata_blk2_results1 %>%
    dplyr::select( Xblk:LATE.hat, SE.wald, n, n_comp, k )


# Calculate reduction in SE for different estimators and reduce to LATEs only
iv_key1 <- strata_blk2_results1[18:24,]  %>%
    dplyr::select( Xblk, pi.hat, LATE.hat, SE.wald, n )

unstrat_wald <- strata_blk2_results1[strata_blk2_results1$Xblk=="UNSTRAT", "SE.wald"]

# Reformat results to include test statistics, p-value, and SE reduction
iv_key1 <- mutate( iv_key1,
                   t = LATE.hat / SE.wald,
                   pvalue = 2 * pmin( 1 - pnorm(t), pnorm( t )),
                   percent_red = 100*SE.wald/unstrat_wald)
iv_key1

# Format into table
print(xtable(iv_key1, digits = c(0,0,3,3,4,0, 0,4, 1 )), include.rownames=FALSE )


#Look at information on stratum level results
min(strata_blk2_results1$pi.hat[1:17])
max(strata_blk2_results1$pi.hat[1:17])
min(strata_blk2_results1$LATE.hat[1:17])
max(strata_blk2_results1$LATE.hat[1:17])



##### MCAFE estimation #####

# Following Blackwell & Pashley (2023) we check factorial type MCAFE
# estimates with phone and mailers

cat( "Analyzing MCAFE\n" )

## Functions in script below implement the same IV estimation strategies but
## applied to MCAFE estimates instead of standard binary IV estimates (that ignore the other factors)
source( here::here( "mcafe post-strat.R") )

# Subset to data without missing assignments on phonecalls and mailers
hh_dat2 <- hh_dat2[complete.cases(hh_dat2[,c("phonegrp", "mailgrp")]),]

# Perform analysis using MCAFE
strata_blk2_results2<- IV.est.strat.fac( data = hh_dat2,  Yobs="any_to", S="ipc", Z="persngrp", W="phonegrp", W2 = "mailgrp", strat_var ="blk2" )
tail(strata_blk2_results2)

# Reduce to key results
iv_key2 <- strata_blk2_results2[18:22,]  %>%
    dplyr::select( Xblk, pi.hat, LATE.hat, SE.wald, n )

# Reformat results to include test statistics, p-value, and SE reduction
unstrat_wald2 <- iv_key2[iv_key2$Xblk=="UNSTRAT", "SE.wald"]

iv_key2 <- mutate( iv_key2,
                   t = LATE.hat / SE.wald,
                   pvalue = 2 * pmin( 1 - pnorm(t), pnorm( t )),
                   percent_red = 100*SE.wald/unstrat_wald2 )
iv_key2



