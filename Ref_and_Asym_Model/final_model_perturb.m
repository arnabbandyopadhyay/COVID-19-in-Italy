
clear all
%% Initial population
P0.abruzzo=1311580;P0.campania=5801692;P0.emilia=4459477;P0.friuli=1215220;
P0.lazio=5879082;P0.liguria=1550640;P0.lombardia=10060574;
P0.marche=1525271;P0.molise=305617;P0.piemonte=4356406;P0.sardegna=1639591;
P0.sicilia=4999891;P0.toscana=3729641;P0.valledaosta=125666;P0.veneto=4905854;
P0.calabria=1947131;P0.puglia=4029053;P0.umbria=882015;P0.basilicata=562869;
P0.italy=60359546;P0.bolzano=520891;P0.trento=538223;

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
% comment in case of asymptomatic model

C0.abruzzo=1/(1-0.895);C0.campania=3/(1-0.5);C0.emilia=18/(1-0.9169);C0.friuli=6/(1-0.8902);
C0.italy=221/(1-0.9377);C0.lazio=3/(1-0.5);C0.liguria=1/(1-0.9302);C0.lombardia=225/(1-0.9558);
C0.marche=1/(1-0.9057);C0.molise=3/(1-0.5124);C0.piemonte=2/(1-0.929);C0.sardegna=1/(1-0.9553);
C0.sicilia=1/(1-0.5);C0.toscana=2/(1-0.8777);C0.valledaosta=2/(1-0.9025);C0.veneto=32/(1-0.8967);
C0.calabria=1/(1-0.8416);C0.puglia=1/(1-0.9437);C0.umbria=2/(1-0.5);C0.basilicata=1/(1-0.574);
C0.bolzano=1/(1-0.9256);C0.trento=4/(1-0.9161);

%% Asymptomatic model: \bar{mu}=alpha (mu=1)
% uncomment in case of asymptomatic model
% C0.abruzzo=1/(1-0.4);C0.campania=3/(1-0.4);C0.emilia=18/(1-0.4);C0.friuli=6/(1-0.4);
% C0.italy=221/(1-0.4);C0.lazio=3/(1-0.4);C0.liguria=1/(1-0.4);C0.lombardia=225/(1-0.4);
% C0.marche=1/(1-0.4);C0.molise=3/(1-0.4);C0.piemonte=2/(1-0.4);C0.sardegna=1/(1-0.4);
% C0.sicilia=1/(1-0.4);C0.toscana=2/(1-0.4);C0.valledaosta=2/(1-0.4);C0.veneto=32/(1-0.4);
% C0.calabria=1/(1-0.4);C0.puglia=1/(1-0.4);C0.umbria=2/(1-0.4);C0.basilicata=1/(1-0.4);
% C0.bolzano=1/(1-0.4);C0.trento=4/(1-0.4);

%% areas names

CN.abruzzo='Abruzzo';CN.campania='Campania';CN.emilia='Emilia';CN.friuli='Friuli';
CN.italy='Italy';CN.lazio='Lazio';CN.liguria='Liguria';CN.lombardia='Lombardia';
CN.marche='Marche';CN.molise='Molise';CN.piemonte='Piemonte';CN.sardegna='Sardegna';
CN.sicilia='Sicilia';CN.toscana='Toscana';CN.valledaosta='Valledaosta';CN.veneto='Veneto';
CN.calabria='Calabria';CN.puglia='Puglia';CN.umbria='Umbria';CN.basilicata='Basilicata';
CN.bolzano='Bolzano';CN.trento='Trento';

%% missing data
% Nation wide data is available from 24/02/2020. 
% some regions started counting later 
mis.abruzzo=3; mis.campania=3;mis.emilia=0;mis.friuli=6;
mis.italy=0;mis.lazio=0;mis.liguria=1;mis.lombardia=0;
mis.marche=2;mis.molise=8;mis.piemonte=0;mis.sardegna=8;
mis.sicilia=1;mis.toscana=1;mis.valledaosta=10;mis.veneto=0;
mis.calabria=4;mis.puglia=3;mis.umbria=6;mis.basilicata=8;
mis.bolzano=1;mis.trento=8;

%% Undetected amount
% For most of the regions we estimated the undetected amount through
% Bayesian MCMC framework. For some regions data were not available, so we
% used undetected=0.5 
undet.abruzzo=0.895; undet.campania=0.5; undet.emilia=0.9169; undet.friuli=0.8902;
undet.italy=0.9377; undet.lazio=0.5; undet.liguria=0.9302; undet.lombardia=0.9558;
undet.marche=0.9057; undet.molise=0.5124;undet.piemonte=0.929; undet.sardegna=0.9553;
undet.sicilia=0.5; undet.toscana=0.8777; undet.valledaosta=0.9025;undet.veneto=0.8967;
undet.calabria=0.8416;undet.puglia=0.9437;undet.umbria=0.5;undet.basilicata=0.574;
undet.bolzano=0.9256; undet.trento=0.9161;
 
%%
% toevaluate=fieldnames(P0);
toevaluate={'lombardia'};


for cfile=1:length(toevaluate)

    copyfile('Models/last_model_t_undetected.def',...
        sprintf('Models/last_model_t_undetected_%s.def',getfield(CN,toevaluate{cfile})));
    copyfile('Models/last_model_t_undetected_loop.def',...
        sprintf('Models/last_model_t_undetected_loop_%s.def',getfield(CN,toevaluate{cfile})));


end




parfor key=1:length(toevaluate)
    
     disp(key);
    
     %% POI library path: need to update the path accordingly
     
     javaaddpath('/Users/Arnab/projects/poi_library/poi-3.8-20120326.jar');
     javaaddpath('/Users/Arnab/projects/poi_library/poi-ooxml-3.8-20120326.jar');
     javaaddpath('/Users/Arnab/projects/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
     javaaddpath('/Users/Arnab/projects/poi_library/xmlbeans-2.3.0.jar');
     javaaddpath('/Users/Arnab/projects/covid-19/poi_library/dom4j-1.6.1.jar');
     

     

     %% declare the home directory: needs to update accordingly
     
    cd '/Users/Arnab/Desktop/public_code/Ref_and_Asym_Model'
    
    
    %% modifying Susceptible and Exposed equations and initial conditions
    
    line1=sprintf('"-r1*(x3+x4+ki*x14+bet*(x5+x6))*x1/%d"     // sus x1',getfield(P0,toevaluate{key}));
    line2=sprintf('"r1*(x3+x4+ki*x14+bet*(x5+x6))*x1/%d -(r3/(r3*5.2 -1))*x2"     // exposed x2',getfield(P0,toevaluate{key}));
    init_x1=sprintf('init_x1      "%d"',getfield(P0,toevaluate{key}));
    init_x2=sprintf('init_x2      "%d"',getfield(C0,toevaluate{key}));
    
    Func_replace_string(sprintf('Models/last_model_t_undetected_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line1',line1);
   
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'line2',line2);
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'init_x1',init_x1);
    
    Func_replace_string(sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/covid_19_main_2_%s.def',getfield(CN,toevaluate{key})),'init_x2',init_x2);
    
    %loop file modification
    
    Func_replace_string(sprintf('Models/last_model_t_undetected_loop_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line1',line1);
    
    Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key})),'line2',line2);

    
    %% Load the data
    
    datatoread=sprintf('Data/%s.xls',getfield(CN,toevaluate{key}));

    %% Load D2D framework  
    
%     para=load(sprintf('/home/msi18/Desktop/covid-19/delta_folder/NOrec/area_specific_%s.mat',getfield(CN,toevaluate{key})));
    
    %% Reference model
    undetected=getfield(undet,toevaluate{key});
    
    %% Asymptomatic model
    %uncomment if asymptomatic model is running: 
    %overwrites undetected amount
%     undetected=0.4;
    
    %% load model file and window model file
    
    deffile=sprintf('covid_19_main_2_%s',getfield(CN,toevaluate{key}));
    loopfile=sprintf('last_model_t_loop_2_%s.def',getfield(CN,toevaluate{key}));
    
    
    cityname=getfield(CN,toevaluate{key});
    missing=getfield(mis,toevaluate{key});
%     ar=firstloop(datatoread,deffile,loopfile,cityname,...
%         missing, undetected);%,para.savpars);
    ar=perturb_loop(datatoread,deffile,loopfile,cityname,...
        missing, undetected);%,para.savpars);

    
    delete *mexmaci64



end

