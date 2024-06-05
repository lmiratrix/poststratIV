

library( poststratIV )

test_that("Testing the simulation functions", {



        # Making data works?
        make_dat( 2000, 0, 0.5, pi_c = 0.05,
                  pred_comp = "no", pred_Y = "no", het_tx="yes",
                  params = default_sim_params() )

        set.seed( 40404 )
        a <- poststratIV:::one_run( 2000, 0, 0.5, pi_c = 0.05,
                      pred_comp = "no", pred_Y = "no", het_tx="yes" )

        # Should be the same as above
        set.seed( 40404 )
        b <- poststratIV:::one_run( 2000, 0, 0.5, pi_c = 0.05,
                      pred_comp = "no", pred_Y = "no", het_tx="yes",
                      params = default_sim_params() )

        expect_equal( a, b )


        rs <- run_sim( reps = 4, N = 2000, pi_c = 0.05,
                 nt_shift = -1, sd0 = 1.0, pred_comp = "no",
                 pred_Y = "no", het_tx = "yes" )

        rs

        #rs2 <- run_sim( reps = 1000, N = 500, pi_c = 0.05, nt_shift = -1, het_tx = "no",
        #         sd0 = 1.0, pred_comp = "no", pred_Y = "no" )

        #rs2

        # two-sided checking
        set.seed( 40440 )
        ts_dat <- make_dat( 2000, pi_c = 0.05, pred_comp = "no", pred_Y = "no",
                  params = default_sim_params(), one_sided = FALSE )

        head( ts_dat )
        tb = table( ts_dat$S0, ts_dat$S1 )
        tb
        expect_true( sum( tb != 0 ) == 3 )
})
