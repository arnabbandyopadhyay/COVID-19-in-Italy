# Reference and Asymptomatic models

**Steps to generate Reference Model single simulation:**

After proper setting up POI library in MATLAB:
1) run 'setup_reference.m' file
2) run the following:
  mkdir Ref_results
  mv *.mat Ref_results
  mv *.csv Ref_results
  mv Ref_results Perturbation

**Steps to generate Asymptomatic Model single simulation:**
After proper setting up POI library in MATLAB:
1) run 'setup_asymptomatic.m' file
2) run the following:
  mkdir Asymp_results
  mv *.mat Asymp_results
  mv *.csv Asymp_results
  mv Asymp_results Perturbation

**Steps for perturbation analysis:**
1) run the following:

    cd Perturbation
    
    run 'ReferPert.m' --> it generates perturbed simulation for Reference and Testing(Early&Late) Model
    
    run 'AsympPert.m' --> it generates perturbed simulations for Asymptomatic Model

**Steps to generate plots:**
1) run 'DOPAR_fit_rt.R' --> it generates main and supplementary fit(from Ref.Mod.) and Rt plots(Ref and Asymp comparison)

2) run 'DOPAR_Testing.R' --> it generates main and supplementary figures for the Early and Late Tesing Model

