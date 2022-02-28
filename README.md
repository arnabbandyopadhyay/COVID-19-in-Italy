# Testing and isolation to prevent overloaded health care facilities and to reduce death rates in the SARS-CoV-2 pandemic in Italy

Arnab Bandyopadhyay*, Marta Schips*, Tanmay Mitra, Sahamoddin Khailaie, Sebastian Binder, Michael Meyer-Hermann*

Department of Systems Immunology and Braunschweig Integrated Centre of Systems Biology(BRICS), Helmholtz Centre for Infection Research, Braunschweig, Germany

*Correspondence:
mmh@theoretical-biology.de 
arnab.bandyopadhyay@theoretical-biology.de 
marta.schips@theoretical-biology.de

https://www.medrxiv.org/content/10.1101/2020.10.12.20211169v1

 **system requirements:**
 
MATLAB (best R2018b or newer)

MATLAB Symbolic Toolbox

MATLAB Optimization Toolbox

MATLAB Parallel Computing Toolbox 

Data2Dynamics (add the folder /arFramework3 to the MATLAB path). For installation check: 
https://github.com/Data2Dynamics/d2d/wiki/Installation

POI library

**Note:**
It takes approximately one one hour (on Intel Xeon processor @ 3.30GHz) for a single region and for the full dataset. The code is already parallelized in MATLAB and therefore depending upon the availability of cores multiple regions can be evaluated simultaneously.

**Folder organization**

PARALLEL folder contains two folders, Capacity and Reference_Asymptomatic. Reference_Asymptomatic folder contains all necessary files required to simulate Reference model, Asymptomatic model and Testing model (early or late). Data folder contains regional and national data on active cases (‘qua’), hospitalized (‘hos’), ICU (‘icu’) and death (‘dead’) numbers to calibrate the model. Models folder contains model files necessary to evaluate the model in Data2Dynamics framework. Perturbation folder contains main results and matlab scripts to generate perturbation. Ref_results and Asymp_results contain parameter files for each region generated by Reference model and Asymptomatic model, respectively. Refer_pert_res, Asym_pert_res_testing_pert_res contains perturbation result of Reference, Asymptomatic and Testing model respectively.  Refer_rt_res, Asym_rt_res_testing_rt_res contain perturb Rt values of Reference, Asymptomatic and Testing model respectively. 

**Important to keep in mind**

a) Model equations, observables and initial conditions are declared in the model files in Models folder. For example, 'model_Italy.def' is a template file where susceptible and exposed state equations and corresponding initial conditions are empty as the susceptible and exposed amount is different for different regions and are filled ('model_2_Italy.def') while evaluating the relevant .m file. 'model_2_*region_name*.def' was used to estimate the parameters by fitting the first 14 days of data points considering exponential growth. 

‘model_loop_2_*region_name*.def’ is a template file where all initial conditions are empty as the initial conditions of current window is filled by the state space of the previous window and is modified on the way. This file was used in windowing procedure.

'setup_reference.m' file in Reference_Asymptomatic folder will modify the template files and continue parameter estimation process by calling 'firstloop.m' function.

Similarly, 'setup_asymptomatic.m' file in Reference_Asymptomatic folder can be used to evaluate Asymptomatic model.

Perturbations were performed in Perturbation folder. Additional instructions are given in the Reference_Asymptomatic and Capacity folder on how to generate results, moving it in right directory and running the Rscript to generate figures 4-7 in the main text and supplementary figures 7-13. 

b) POI library path needs to be updated accordingly. 

c) home directory must contain following functions: ‘R0calc.m’, ‘xlwrite.m’ and ‘Func_replace_string.m’

Additional instructions are given in the readme file in the folder to generate figures

Data folder contain data of all regions as of 23 July, 2020 and downloaded from https://github.com/pcm-dpc/COVID-19

