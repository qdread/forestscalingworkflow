# All CmdStan calls for fitting piecewise models (on slurm cluster)
# QDR / Forestlight / 29 October 2019

# All fits allow the user to specify random seed and number of iterations.

# All scalings: DIAMETER (x) AND BIOMASS PRODUCTION (y)
# =====================================================

# 1 segment density

sbatch --export=scaling=production,guild=fg1,year=1995,model=density1,NS=1000,NW=5000,seed=11 --time=4:00:00 --job-name=d1_1 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg2,year=1995,model=density1,NS=1000,NW=5000,seed=12 --time=4:00:00 --job-name=d1_2 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg3,year=1995,model=density1,NS=1000,NW=5000,seed=13 --time=4:00:00 --job-name=d1_3 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg4,year=1995,model=density1,NS=1000,NW=5000,seed=14 --time=4:00:00 --job-name=d1_4 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg5,year=1995,model=density1,NS=1000,NW=5000,seed=15 --time=4:00:00 --job-name=d1_5 fitpiecewise.sh
sbatch --export=scaling=production,guild=unclassified,year=1995,model=density1,NS=1000,NW=5000,seed=16 --time=4:00:00 --job-name=d1_u fitpiecewise.sh
sbatch --export=scaling=production,guild=alltree,year=1995,model=density1,NS=1000,NW=5000,seed=17 --time=4:00:00 --job-name=d1_a fitpiecewise.sh

# 2 segment density

sbatch --export=scaling=production,guild=fg1,year=1995,model=density2,NS=1000,NW=5000,seed=211 --time=4:00:00 --job-name=d2_1 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg2,year=1995,model=density2,NS=1000,NW=5000,seed=212 --time=4:00:00 --job-name=d2_2 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg3,year=1995,model=density2,NS=1000,NW=5000,seed=213 --time=7-00:00:00 --job-name=d2_3 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg4,year=1995,model=density2,NS=1000,NW=5000,seed=214 --time=4:00:00 --job-name=d2_4 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg5,year=1995,model=density2,NS=1000,NW=5000,seed=215 --time=4:00:00 --job-name=d2_5 fitpiecewise.sh
sbatch --export=scaling=production,guild=unclassified,year=1995,model=density2,NS=1000,NW=5000,seed=216 --time=4:00:00 --job-name=d2_u fitpiecewise.sh
sbatch --export=scaling=production,guild=alltree,year=1995,model=density2,NS=1000,NW=5000,seed=217 --time=7-00:00:00 --job-name=d2_a fitpiecewise.sh

# 3 segment density

sbatch --export=scaling=production,guild=fg1,year=1995,model=density3,NS=1000,NW=5000,seed=311 --time=4:00:00 --job-name=d3_1 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg2,year=1995,model=density3,NS=1000,NW=5000,seed=312 --time=4:00:00 --job-name=d3_2 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg3,year=1995,model=density3,NS=1000,NW=5000,seed=313 --time=7-00:00:00 --job-name=d3_3 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg4,year=1995,model=density3,NS=1000,NW=5000,seed=314 --time=4:00:00 --job-name=d3_4 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg5,year=1995,model=density3,NS=1000,NW=5000,seed=315 --time=4:00:00 --job-name=d3_5 fitpiecewise.sh
sbatch --export=scaling=production,guild=unclassified,year=1995,model=density3,NS=1000,NW=5000,seed=316 --time=4:00:00 --job-name=d3_u fitpiecewise.sh
sbatch --export=scaling=production,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=317 --time=7-00:00:00 --job-name=d3_a fitpiecewise.sh

# 1 segment production

sbatch --export=scaling=production,guild=fg1,year=1995,model=production1,NS=1000,NW=5000,seed=411 --time=4:00:00 --job-name=p1_1 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg2,year=1995,model=production1,NS=1000,NW=5000,seed=412 --time=4:00:00 --job-name=p1_2 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg3,year=1995,model=production1,NS=1000,NW=5000,seed=413 --time=4:00:00 --job-name=p1_3 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg4,year=1995,model=production1,NS=1000,NW=5000,seed=414 --time=4:00:00 --job-name=p1_4 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg5,year=1995,model=production1,NS=1000,NW=5000,seed=415 --time=4:00:00 --job-name=p1_5 fitpiecewise.sh
sbatch --export=scaling=production,guild=unclassified,year=1995,model=production1,NS=1000,NW=5000,seed=416 --time=4:00:00 --job-name=p1_u fitpiecewise.sh
sbatch --export=scaling=production,guild=alltree,year=1995,model=production1,NS=1000,NW=5000,seed=417 --time=4:00:00 --job-name=p1_a fitpiecewise.sh

# 2 segment production

sbatch --export=scaling=production,guild=fg1,year=1995,model=production2,NS=1000,NW=5000,seed=511 --time=4:00:00 --job-name=p2_1 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg2,year=1995,model=production2,NS=1000,NW=5000,seed=512 --time=4:00:00 --job-name=p2_2 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg3,year=1995,model=production2,NS=1000,NW=5000,seed=513 --time=7-00:00:00 --job-name=p2_3 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg4,year=1995,model=production2,NS=1000,NW=5000,seed=514 --time=7-00:00:00 --job-name=p2_4 fitpiecewise.sh
sbatch --export=scaling=production,guild=fg5,year=1995,model=production2,NS=1000,NW=5000,seed=515 --time=4:00:00 --job-name=p2_5 fitpiecewise.sh
sbatch --export=scaling=production,guild=unclassified,year=1995,model=production2,NS=1000,NW=5000,seed=516 --time=4:00:00 --job-name=p2_u fitpiecewise.sh
sbatch --export=scaling=production,guild=alltree,year=1995,model=production2,NS=1000,NW=5000,seed=517 --time=7-00:00:00 --job-name=p2_a fitpiecewise.sh

# Individual power laws: DIAMETER (x) AND RAW INCOMING LIGHT ENERGY (y)
# =====================================================================

# 1 segment

sbatch --export=scaling=rawlightscaling,guild=fg1,year=1995,model=production1,NS=1000,NW=5000,seed=1011 --time=4:00:00 --job-name=l1_1 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg2,year=1995,model=production1,NS=1000,NW=5000,seed=1012 --time=4:00:00 --job-name=l1_2 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg3,year=1995,model=production1,NS=1000,NW=5000,seed=1013 --time=4:00:00 --job-name=l1_3 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg4,year=1995,model=production1,NS=1000,NW=5000,seed=1014 --time=4:00:00 --job-name=l1_4 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg5,year=1995,model=production1,NS=1000,NW=5000,seed=1015 --time=4:00:00 --job-name=l1_5 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=unclassified,year=1995,model=production1,NS=1000,NW=5000,seed=1016 --time=4:00:00 --job-name=l1_u fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=alltree,year=1995,model=production1,NS=1000,NW=5000,seed=1017 --time=4:00:00 --job-name=l1_a fitpiecewise.sh

# 2 segment

sbatch --export=scaling=rawlightscaling,guild=fg1,year=1995,model=production2,NS=1000,NW=5000,seed=1111 --time=4:00:00 --job-name=l2_1 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg2,year=1995,model=production2,NS=1000,NW=5000,seed=1112 --time=4:00:00 --job-name=l2_2 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg3,year=1995,model=production2,NS=1000,NW=5000,seed=1113 --time=7-00:00:00 --job-name=l2_3 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg4,year=1995,model=production2,NS=1000,NW=5000,seed=1114 --time=7-00:00:00 --job-name=l2_4 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=fg5,year=1995,model=production2,NS=1000,NW=5000,seed=1115 --time=4:00:00 --job-name=l2_5 fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=unclassified,year=1995,model=production2,NS=1000,NW=5000,seed=1116 --time=4:00:00 --job-name=l2_u fitpiecewise.sh
sbatch --export=scaling=rawlightscaling,guild=alltree,year=1995,model=production2,NS=1000,NW=5000,seed=1117 --time=7-00:00:00 --job-name=l2_a fitpiecewise.sh

# Individual power laws: DIAMETER (x) AND DIAMETER GROWTH RATE (y)
# ================================================================

# 1 segment

sbatch --export=scaling=diamgrowthscaling,guild=fg1,year=1995,model=production1,NS=1000,NW=5000,seed=6611 --time=4:00:00 --job-name=dg1_1 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg2,year=1995,model=production1,NS=1000,NW=5000,seed=6612 --time=4:00:00 --job-name=dg1_2 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg3,year=1995,model=production1,NS=1000,NW=5000,seed=6613 --time=7-00:00:00 --job-name=dg1_3 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg4,year=1995,model=production1,NS=1000,NW=5000,seed=6614 --time=7-00:00:00 --job-name=dg1_4 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg5,year=1995,model=production1,NS=1000,NW=5000,seed=6615 --time=4:00:00 --job-name=dg1_5 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=unclassified,year=1995,model=production1,NS=1000,NW=5000,seed=6616 --time=4:00:00 --job-name=dg1_u fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=alltree,year=1995,model=production1,NS=1000,NW=5000,seed=6617 --time=7-00:00:00 --job-name=dg1_a fitpiecewise.sh

# 2 segment

sbatch --export=scaling=diamgrowthscaling,guild=fg1,year=1995,model=production2,NS=1000,NW=5000,seed=5511 --time=4:00:00 --job-name=dg2_1 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg2,year=1995,model=production2,NS=1000,NW=5000,seed=5512 --time=4:00:00 --job-name=dg2_2 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg3,year=1995,model=production2,NS=1000,NW=5000,seed=5513 --time=7-00:00:00 --job-name=dg2_3 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg4,year=1995,model=production2,NS=1000,NW=5000,seed=5514 --time=7-00:00:00 --job-name=dg2_4 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=fg5,year=1995,model=production2,NS=1000,NW=5000,seed=5515 --time=4:00:00 --job-name=dg2_5 fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=unclassified,year=1995,model=production2,NS=1000,NW=5000,seed=5516 --time=4:00:00 --job-name=dg2_u fitpiecewise.sh
sbatch --export=scaling=diamgrowthscaling,guild=alltree,year=1995,model=production2,NS=1000,NW=5000,seed=5517 --time=7-00:00:00 --job-name=dg2_a fitpiecewise.sh

# Individual power laws: DIAMETER (x) and CROWN VOLUME (y)
# ========================================================

# 1 segment

sbatch --export=scaling=volumescaling,guild=fg1,year=1995,model=production1,NS=1000,NW=5000,seed=9911 --time=4:00:00 --job-name=v1_1 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg2,year=1995,model=production1,NS=1000,NW=5000,seed=9912 --time=4:00:00 --job-name=v1_2 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg3,year=1995,model=production1,NS=1000,NW=5000,seed=9913 --time=7-00:00:00 --job-name=v1_3 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg4,year=1995,model=production1,NS=1000,NW=5000,seed=9914 --time=7-00:00:00 --job-name=v1_4 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg5,year=1995,model=production1,NS=1000,NW=5000,seed=9915 --time=4:00:00 --job-name=v1_5 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=unclassified,year=1995,model=production1,NS=1000,NW=5000,seed=9916 --time=7-00:00:00 --job-name=v1_u fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=alltree,year=1995,model=production1,NS=1000,NW=5000,seed=9917 --time=7-00:00:00 --job-name=v1_a fitpiecewise.sh

# 2 segment

sbatch --export=scaling=volumescaling,guild=fg1,year=1995,model=production2,NS=1000,NW=5000,seed=1911 --time=4:00:00 --job-name=v2_1 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg2,year=1995,model=production2,NS=1000,NW=5000,seed=1912 --time=4:00:00 --job-name=v2_2 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg3,year=1995,model=production2,NS=1000,NW=5000,seed=1913 --time=7-00:00:00 --job-name=v2_3 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg4,year=1995,model=production2,NS=1000,NW=5000,seed=1914 --time=7-00:00:00 --job-name=v2_4 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=fg5,year=1995,model=production2,NS=1000,NW=5000,seed=1915 --time=4:00:00 --job-name=v2_5 fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=unclassified,year=1995,model=production2,NS=1000,NW=5000,seed=1916 --time=7-00:00:00 --job-name=v2_u fitpiecewise.sh
sbatch --export=scaling=volumescaling,guild=alltree,year=1995,model=production2,NS=1000,NW=5000,seed=1917 --time=7-00:00:00 --job-name=v2_a fitpiecewise.sh

# All groups together in multilevel model: LIGHT PER AREA (x) and MORTALITY (y)
# =============================================================================

sbatch --job-name=mort fitmortality.sh

# Von Bertalanffy model: LIGHT PER AREA (x) and PRODUCTION PER AREA (y)
# =====================================================================

sbatch --export=scaling=light,guild=fg1,year=1995,model=vonb,NS=1000,NW=5000,seed=8811 --time=4:00:00 --job-name=vonb_1 fitpiecewise.sh
sbatch --export=scaling=light,guild=fg2,year=1995,model=vonb,NS=1000,NW=5000,seed=8812 --time=4:00:00 --job-name=vonb_2 fitpiecewise.sh
sbatch --export=scaling=light,guild=fg3,year=1995,model=vonb,NS=1000,NW=5000,seed=8813 --time=4:00:00 --job-name=vonb_3 fitpiecewise.sh
sbatch --export=scaling=light,guild=fg4,year=1995,model=vonb,NS=1000,NW=5000,seed=8814 --time=4:00:00 --job-name=vonb_4 fitpiecewise.sh
sbatch --export=scaling=light,guild=fg5,year=1995,model=vonb,NS=1000,NW=5000,seed=8815 --time=4:00:00 --job-name=vonb_5 fitpiecewise.sh
sbatch --export=scaling=light,guild=unclassified,year=1995,model=vonb,NS=1000,NW=5000,seed=8816 --time=4:00:00 --job-name=vonb_u fitpiecewise.sh
sbatch --export=scaling=light,guild=alltree,year=1995,model=vonb,NS=1000,NW=5000,seed=8817 --time=4:00:00 --job-name=vonb_a fitpiecewise.sh

# Chains that did not converge (rerun)
# ====================================

sbatch --export=scaling=production,guild=alltree,year=1995,model=production2,NS=1000,NW=5000,seed=2517 --time=7-00:00:00 --job-name=p2_a --array=2 fitpiecewise.sh

sbatch --export=scaling=production,guild=fg3,year=1995,model=density3,NS=1000,NW=5000,seed=3 --time=4-00:00:00 --job-name=d3_3 fitpiecewise_newargs.sh # old model, new model fitting parameters

sbatch --export=scaling=production,guild=fg3,year=1995,model=density3_decreasingslopes,NS=1000,NW=7500,seed=666 --time=4-00:00:00 --job-name=d3_newprior fitpiecewise.sh # model with tightened priors.
sbatch --export=scaling=production,guild=fg3,year=1995,model=density3_decreasingslopes,NS=1000,NW=7500,seed=888 --time=4-00:00:00 --job-name=d3_newprior --array=4-30 fitpiecewise.sh # model with tightened priors. # Use 8,21,25.
#sbatch --export=scaling=production,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=11317 --time=7-00:00:00 --job-name=d3_a --array=1-2 fitpiecewise.sh
#sbatch --export=scaling=volumescaling,guild=unclassified,year=1995,model=production1,NS=1000,NW=5000,seed=19916 --time=7-00:00:00 --job-name=v1_u fitpiecewise.sh
#sbatch --export=scaling=diamgrowthscaling,guild=fg4,year=1995,model=production2,NS=1000,NW=5000,seed=9514 --time=7-00:00:00 --job-name=dg2_4 --array=2 fitpiecewise.sh

#sbatch --export=scaling=production,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=33317 --time=7-00:00:00 --job-name=d3_a --array=2 fitpiecewise.sh


