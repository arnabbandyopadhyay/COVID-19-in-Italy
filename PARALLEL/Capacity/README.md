#Capacity models#

--> generate Capacity Model single simulation:
-- run 'exp_cap.m' file
-- run the following:
mkdir CapacResults
mv *.mat CapacResults
mv *.csv CapacResults
mv CapacResults Perturbation

--> perturbation:
run the following:
cd Perturbation
-- run 'captTest_perD.m' --> it generates perturbed simulation for Capacity, MaxCap and Testing(Early&Late) Models

--> plots:
-- run 'DOPAR_CapDeCap.R' --> it generates main and supplementary fit of the Capacity Model and the Excess Dead
-- run 'DOPAR_Testing02.R' --> it generates main and supplementary figures for the Early and Late Tesing Model
-- run 'DOPAR_Testing005.R' --> it generates main and supplementary figures for the Early and Late Tesing Model


TO GENERATE ALPHA PERTURBATION:
-- run 'exp_cap.m' file with different alpha values
for each run with each alpha value
-- run the following:
mkdir alpha*
mv *.mat alpha*
mv *.csv alpha
mv alpha Perturbation
e.g. alpha=0.2 --> *=02

--> perturbation:
run the following:
cd Perturbation
-- run 'AlphaPerturb.m' --> WARNING: change the ALPHA variable at the beginning of the script with the name of the generated alpha* folders
-- run 'onlyDead_15122021.R' --> WARNING: change the alphaAA variable at the beginning of the script with the name of the generate alpha* folders
