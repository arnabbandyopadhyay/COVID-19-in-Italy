% javaaddpath('/home/msi18/Desktop/covid-19/poi_library/poi-3.8-20120326.jar');
% javaaddpath('/home/msi18/Desktop/covid-19/poi_library/poi-ooxml-3.8-20120326.jar');
% javaaddpath('/home/msi18/Desktop/covid-19/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
% javaaddpath('/home/msi18/Desktop/covid-19/poi_library/xmlbeans-2.3.0.jar');
% javaaddpath('/home/msi18/Desktop/covid-19/poi_library/dom4j-1.6.1.jar');

clear all
home_dir=pwd;

%% set up P0


P0.abruzzo=1311580;P0.campania=5801692;P0.emilia=4459477;P0.friuli=1215220;
P0.lazio=5879082;P0.liguria=1550640;P0.lombardia=10060574;P0.italy=60359546;
P0.marche=1525271;P0.molise=305617;P0.piemonte=4356406;P0.sardegna=1639591;
P0.sicilia=4999891;P0.toscana=3729641;P0.valledaosta=125666;P0.veneto=4905854;
P0.calabria=1947131;P0.puglia=4029053;P0.umbria=882015;P0.basilicata=562869;
P0.bolzano=520891;P0.trento=538223;


%%%initial conditions: undetected model
C0.abruzzo=1*1/(1-0.8939);C0.campania=3*1/(1-0.1856);C0.emilia=18*1/(1-0.9169);C0.friuli=6*1/(1-0.891);
C0.italy=221*1/(1-0.9357);C0.lazio=3*1/(1-0.5);C0.liguria=1*1/(1-0.9309);C0.lombardia=225*1/(1-0.9553);
C0.marche=1*1/(1-0.9057);C0.molise=3*1/(1-0.5186);C0.piemonte=2*1/(1-0.9294);C0.sardegna=1*1/(1-0.9553);
C0.sicilia=1*1/(1-0.5);C0.toscana=2*1/(1-0.8777);C0.valledaosta=2*1/(1-0.9025);C0.veneto=32*1/(1-0.8967);
C0.calabria=1*1/(1-0.8442);C0.puglia=1*1/(1-0.9433);C0.umbria=2*1/(1-0.5);C0.basilicata=1*1/(1-0.5733);
C0.bolzano=1*1/(1-0.9261);C0.trento=4*1/(1-0.9161);

CN.abruzzo='Abruzzo';CN.campania='Campania';CN.emilia='Emilia';CN.friuli='Friuli';
CN.italy='Italy';CN.lazio='Lazio';CN.liguria='Liguria';CN.lombardia='Lombardia';
CN.marche='Marche';CN.molise='Molise';CN.piemonte='Piemonte';CN.sardegna='Sardegna';
CN.sicilia='Sicilia';CN.toscana='Toscana';CN.valledaosta='Valledaosta';CN.veneto='Veneto';
CN.calabria='Calabria';CN.puglia='Puglia';CN.umbria='Umbria';CN.basilicata='Basilicata';
CN.bolzano='Bolzano';CN.trento='Trento';

%%%days missing in the dataset
mis.abruzzo=3; mis.campania=3;mis.emilia=0;mis.friuli=6;
mis.italy=0;mis.lazio=0;mis.liguria=1;mis.lombardia=0;
mis.marche=2;mis.molise=8;mis.piemonte=0;mis.sardegna=8;
mis.sicilia=1;mis.toscana=1;mis.valledaosta=10;mis.veneto=0;
mis.calabria=4;mis.puglia=3;mis.umbria=6;mis.basilicata=8;
mis.bolzano=1;mis.trento=8;

%%%undetected % calculated from IFR
undet.abruzzo=0.8939; undet.campania=0.1856; undet.emilia=0.9169; undet.friuli=0.891;
undet.italy=0.9357; undet.lazio=0.5; undet.liguria=0.9309; undet.lombardia=0.9553;
undet.marche=0.9057; undet.molise=0.5186;undet.piemonte=0.9294; undet.sardegna=0.9553;
undet.sicilia=0.5; undet.toscana=0.8777; undet.valledaosta=0.9025;undet.veneto=0.8967;
undet.calabria=0.8442;undet.puglia=0.9433;undet.umbria=0.5;undet.basilicata=0.5733;
undet.bolzano=0.9261; undet.trento=0.9161;
%%%
%%
%key=8;
toevaluate=fieldnames(P0);

% toevaluate={'bolzano','calabria','campania','friuli','lazio','lombardia','puglia','trento','sardegna','toscana'};

% toevaluate={'italy','lombardia','veneto'};

parfor cfile=1:length(toevaluate)
    
    copyfile('Models/model.def',...
        sprintf('Models/model_%s.def',getfield(CN,toevaluate{cfile})));
    copyfile('Models/model_loop.def',...
        sprintf('Models/model_loop_%s.def',getfield(CN,toevaluate{cfile})));
    
    
end

parfor key=1:length(toevaluate)
    disp(key);
    
    javaaddpath('/home/abp19/Projects/Covid/parallel_it/poi_library/poi-3.8-20120326.jar');
    javaaddpath('/home/abp19/Projects/Covid/parallel_it/poi_library/poi-ooxml-3.8-20120326.jar');
    javaaddpath('/home/abp19/Projects/Covid/parallel_it/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
    javaaddpath('/home/abp19/Projects/Covid/parallel_it/poi_library/xmlbeans-2.3.0.jar');
    javaaddpath('/home/abp19/Projects/Covid/parallel_it/poi_library/dom4j-1.6.1.jar');
    
    
    cd(home_dir) %'/home/abp19/ITALYCOVID15122021/Ref_Asymp/'
    
    
    line1=sprintf('"-r1*(x3+x4+ki*x14+bet*(x5+x6))*x1/%d"     // sus x1',getfield(P0,toevaluate{key}));
    line2=sprintf('"r1*(x3+x4+ki*x14+bet*(x5+x6))*x1/%d -(r3/(r3*5.2 -1))*x2"     // exposed x2',getfield(P0,toevaluate{key}));
    init_x1=sprintf('init_x1      "%d"',getfield(P0,toevaluate{key}));
    init_x2=sprintf('init_x2      "%d"',getfield(C0,toevaluate{key}));
    
    Copy_of_func_replace_string(sprintf('Models/model_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),'line1',line1);
    
    Copy_of_func_replace_string(sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),'line2',line2);
   
    Copy_of_func_replace_string(sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),'init_x1',init_x1);
   
    Copy_of_func_replace_string(sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/model_2_%s.def',getfield(CN,toevaluate{key})),'init_x2',init_x2);
    
    %loop file modification
    
    Copy_of_func_replace_string(sprintf('Models/model_loop_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/model_loop_2_%s.def',getfield(CN,toevaluate{key})),'line1',line1);
    
    Copy_of_func_replace_string(sprintf('Models/model_loop_2_%s.def',getfield(CN,toevaluate{key})),...
        sprintf('Models/model_loop_2_%s.def',getfield(CN,toevaluate{key})),'line2',line2);
    %     datatoread=sprintf('Data/dat_%s_d2d_tanmay.xls',getfield(CN,toevaluate{key}));
    datatoread=sprintf('../Data/%s.xls',getfield(CN,toevaluate{key}));
    
    
    
    %% Load D2D framework
    
    %%% to change in case of asymptomatic model!
    undetected=getfield(undet,toevaluate{key});
    
    
    deffile=sprintf('model_2_%s',getfield(CN,toevaluate{key}));
    loopfile=sprintf('model_loop_2_%s.def',getfield(CN,toevaluate{key}));
    cityname=getfield(CN,toevaluate{key});
    missing=getfield(mis,toevaluate{key});
    ar=firstloop(datatoread,deffile,loopfile,cityname,...
        missing, undetected,home_dir);
    
    
    %% Get necessary things
    
    
    
    delete *mexmaci64
    
    
    
end
