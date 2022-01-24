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

**Important to keep in mind**

a) Folder names are according to the model name and contains Model and Data folder. Models folder contains two model files, ‘last_model_t_undetected.def’ and ‘last_model_t_undetected_loop.def’.
Model equations, observables and initial conditions are declared in the model files. ‘last_model_t_undetected.def’ is used to estimate the parameters by fitting the first 14 days of data points considering exponential growth. Susceptible and exposed state equations and corresponding initial conditions are empty as the susceptible and exposed amount is different for different regions and are filled while evaluating the relevant .m file

‘last_model_t_undetected_loop.def’ is used in windowing and all initial conditions are empty as the initial conditions of current window is filled by the state space of the previous window and is modified on the way 

b) POI library path and home directory needs to be updated accordingly 

c) home directory must contain following functions: ‘R0calc.m’, ‘xlwrite.m’ and ‘Func_replace_string.m’

Additional instructions are given in the readme file in the folder to generate figures

Data_upto_23_July contain data of all regions as of 23 July, 2020 and downloaded from https://github.com/pcm-dpc/COVID-19

