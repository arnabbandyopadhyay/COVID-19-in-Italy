
clear all

%% population
P0.lombardia=10060574;P0.italy=60359546;P0.veneto=4905854;

%% Initial exposed
% Initial exposed persons are the first reported amount of infected people in the data.
% we assumed that these amount of infected people were exposed before the
% incubation period (1/R2 + 1/R3) and we start the simulation from -incubation period so
% that these exposed persons converted into the infection after incubation period. 
% As our model differentiates between detected and undetected, 
% exposed amount was further calibrated by the amount of undetected present 
% (estimated through Bayesian MCMC).
% According to the model, mu * (1-alpha) * Exposed = reported infection (first entry in data)
% Exposed= reported infection/(mu*(1-alpha)); mu*(1-alpha)= 1-\bar{mu}
% which transforms into:  Exposed = reported infection/(1-\bar{mu})

C0.italy=221/(1-0.9377);C0.lombardia=225/(1-0.9558);C0.veneto=32/(1-0.8967);

%% Areas names
CN.italy='Italy';CN.lombardia='Lombardia';CN.veneto='Veneto';

mis.italy=0;mis.lombardia=0;mis.veneto=0;

%% Undetected amount
% estimated through Bayesian MCMC framework.
undet.italy=0.9357; undet.lombardia=0.9553; undet.veneto=0.8967;

toevaluate={'italy','lombardia','veneto'};

for cfile=1:length(toevaluate)
    
    copyfile('Models/last_model_t_undetected.def',...
        sprintf('Models/last_model_t_undetected_%s.def',getfield(CN,toevaluate{cfile})));
    copyfile('Models/last_model_t_undetected_loop.def',...
        sprintf('Models/last_model_t_undetected_loop_%s.def',getfield(CN,toevaluate{cfile})));
    
    
end


parfor key=1:length(toevaluate)
    
    disp(key);
    
    %% POI library path: need to update accordingly
    javaaddpath('/Users/Arnab/projects/poi_library/poi-3.8-20120326.jar');
    javaaddpath('/Users/Arnab/projects/poi_library/poi-ooxml-3.8-20120326.jar');
    javaaddpath('/Users/Arnab/projects/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
    javaaddpath('/Users/Arnab/projects/poi_library/xmlbeans-2.3.0.jar');
    javaaddpath('/Users/Arnab/projects/covid-19/poi_library/dom4j-1.6.1.jar');
    
    %% declare the home directory: needs to update accordingly
    cd '/Users/Arnab/Desktop/public_code/Testing_Model'
    
    %% modifying Susceptible and Exposed equations and initial conditions
    %% introducing IXD compartment
    
    line1=sprintf('"-r1*(x3+x4+ki*x14+bet*(x5+x6+x15))*x1/%d"     // sus x1',getfield(P0,toevaluate{key}));
    line2=sprintf('"r1*(x3+x4+ki*x14+bet*(x5+x6+x15))*x1/%d -(r3/(r3*5.2 -1))*x2"     // exposed x2',getfield(P0,toevaluate{key}));
    line3=sprintf('"((%d - und)/(1-alpha))*r3*x3 -r4*x15"     // Ixd x15',getfield(undet,toevaluate{key}));
    line4=sprintf('"((1-%d)/(1-alpha))*rho*r3*x3-r6*x5" //inf x5 use r6',getfield(undet,toevaluate{key}));
    line5=sprintf('"((1-%d)/(1-alpha))*(1-rho)*r3*x3-r4*x6" //inf2 x6',getfield(undet,toevaluate{key}));
    
    init_x1=sprintf('init_x1      "%d"',getfield(P0,toevaluate{key}));
    init_x2=sprintf('init_x2      "%d"',getfield(C0,toevaluate{key}));
  
    
    Func_replace_string(sprintf('Models/last_model_t_undetected_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line1',line1);
    
   
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line2',line2);
    
    
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line3',line3);
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line4',line4);
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line5',line5);
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'init_x1',init_x1);
    
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'init_x2',init_x2);
    
    %loop file modification
    
    Func_replace_string(sprintf('Models/last_model_t_undetected_loop_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line1',line1);
    
    Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line2',line2);
    
    Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line3',line3);
    
    Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line4',line4);
    
    Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line5',line5);
    
    datatoread=sprintf('Data/%s.xls',getfield(CN,toevaluate{key}));
    
    
    
    %% Load the parameters generated with the Reference model
    % pre-generated files available
    savpars=load(sprintf('ref_results/area_specific_%s.mat',getfield(CN,toevaluate{key})));
    
    undetected=getfield(undet,toevaluate{key});
    %% load model file and window model file
    deffile=sprintf('covid_19_main_2_%s',getfield(CN,toevaluate{key}));
    loopfile=sprintf('last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key}));
    cityname=getfield(CN,toevaluate{key});
    missing=getfield(mis,toevaluate{key});
    
    fixpars=savpars.savpars;
    ar=fixloop_undet(datatoread,deffile,loopfile,cityname,...
        missing, undetected,savpars.savpars);
    
    delete *mexmaci64
    
    
    
end
