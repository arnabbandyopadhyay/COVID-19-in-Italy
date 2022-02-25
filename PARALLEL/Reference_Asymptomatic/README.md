# Reference and Asymptomatic models

** To generate Reference Model single simulation:**
-- run 'setup_reference.m' file
-- run the following:
mkdir Ref_results
mv *.mat Ref_results
mv *.csv Ref_results
mv Ref_results Perturbation

--> generate Asymptomatic Model single simulation:
-- run 'setup_asymptomatic.m' file
-- run the following:
mkdir Asymp_results
mv *.mat Asymp_results
mv *.csv Asymp_results
mv Asymp_results Perturbation

--> perturbation:
run the following:
cd Perturbation
-- run 'ReferPert.m' --> it generates perturbed simulation for Reference and Testing(Early&Late) Model
-- run 'AsympPert.m' --> it generates perturbed simulations for Asymptomatic Model

--> plots:
-- run 'DOPAR_fit_rt.R' --> it generates main and supplementary fit(from Ref.Mod.) and Rt plots(Ref and Asymp comparison)
-- run 'DOPAR_Testing.R' --> it generates main and supplementary figures for the Early and Late Tesing Model

