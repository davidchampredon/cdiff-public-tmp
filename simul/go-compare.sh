#!/bin/sh

rm -rf simul-mc-*.RData

Rscript run-to-compare.R prm_vax_strategy_0.csv prm_vax_p_0.csv > out-run-0.txt 
Rscript run-to-compare.R prm_vax_strategy_frailty40.csv prm_vax_p_0.csv > out-run-f40.txt 
Rscript run-to-compare.R prm_vax_strategy_frailty80.csv prm_vax_p_0.csv > out-run-f80.txt 

# Rscript run-to-compare.R prm_vax_strategy_frailty.csv prm_vax_p_2.csv > out-run-2.txt 
# Rscript run-to-compare.R prm_vax_strategy_frailty.csv prm_vax_p_5.csv > out-run-5.txt 
# Rscript run-to-compare.R prm_vax_strategy_frailty.csv prm_vax_p_8.csv > out-run-8.txt 

echo "   Digesting..."

Rscript compare-sim.R > out-cmp.txt  
Rscript create-main-figures.R > out-main-figs.txt
Rscript simplify-new.R > out-simplify.txt

echo "   Checking empirical vax eff..."

Rscript check-vaxeff.R 1 > check-ve1.out
Rscript check-vaxeff.R 2 > check-ve1.out

echo "   go-compare completed"

