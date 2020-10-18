function ar=fixloop_undet(datatoread,deffile,loopfile,cityname,missing,undetected,fixpars)

shift=5.2; %incubation time
datapoints=7; %window size

parstosave=15; %parameters

accept=0.5; %final portion of undetected

%% set the time at which increased detection starts
day_test=8-1;

%% create the vector for the undet. reduction to use throughout the sim.
if (accept < undetected)
    last=accept/undetected;
    uu=undetected*[repelem(1,day_test) [1:-0.02:last] repelem(last,500)];
else
    uu=undetected*[repelem(1,1000)];
end

%% create a folder with region name, copy all necessary files and change directory. 
% all computation will happen in this folder
mkdir(strcat(cityname))
mkdir(strcat(strcat(cityname),'/Models'))
mkdir(strcat(strcat(cityname),'/Data'))
copyfile(strcat('Models/',strcat(deffile,'.def')), strcat(strcat(cityname),'/Models'))
copyfile(strcat('Models/',loopfile), strcat(strcat(cityname),'/Models'))
copyfile('xlwrite.m',strcat(cityname))
copyfile('Func_replace_string.m',strcat(cityname))
copyfile('R0calc_IXD.m',strcat(cityname))
%take the raw dataset
[num,txt,raw]=xlsread(datatoread);


cd(strcat(cityname))

seq=[3:length(raw)-datapoints+1];
%     seq=[3:5]; % for testing

savinit=zeros(16,length(seq)+1);
savpars=zeros(parstosave,length(seq)+1);
savchi2=zeros(1,length(seq)+1);
savr0=zeros(1,length(seq)+1);
toplot=zeros(1000,6*(length(seq)+1));

savr0per=zeros(100,length(seq)+1);
savchi2per=zeros(100,length(seq)+1);
savfrac=zeros(100,length(seq)+1);

paramsperturb=zeros(100,length(seq)+1);

first_wind=15;
interval=first_wind-missing;
for i=2:interval+1
    raw(i,1)=mat2cell(cell2mat(raw(i,1))+shift-1,1);
end
colel=size(raw,2);
raw(interval+2:end,1:colel)={[]};

%save the new dataset to fit
xlwrite(sprintf('Data/dat_first_%s.xls',cityname), raw);

%% Load D2D framework
%arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])

datafile=sprintf('dat_first_%s',cityname);

arInit;
arLoadModel(deffile);
arLoadData(datafile);

arCompileAll;

%% parameters are fixed as resulting from the Reference model fit
% except for 'und' undetected portion
arSetPars('alpha', log10(fixpars(1,1)),2,1,log10(0.01),log10(0.4));
arSetPars('bet', log10(fixpars(2,1)),2,1,log10(0.05),log10(0.1));
arSetPars('delta', log10(fixpars(3,1)),2,1,log10(0.3),log10(0.9));

arSetPars('ki', log10(fixpars(4,1)),2,1,log10(10^-8),log10(1));

arSetPars('r1', log10(fixpars(5,1)),2,1,log10(0.001),5);
arSetPars('r10', log10(fixpars(6,1)),2,1,log10(0.1),log10(0.9));
arSetPars('r3', log10(fixpars(7,1)),2,1,log10(1/4.2),log10(2/5.2));
arSetPars('r4', log10(fixpars(8,1)),2,1,log10(1/14),log10(1/7));
arSetPars('r5', log10(fixpars(9,1)),2,1,log10(1/16),log10(1/5));
arSetPars('r6', log10(fixpars(10,1)),2,1,log10(1/7),log10(0.9));
arSetPars('r7', log10(fixpars(11,1)),2,1,log10(1/3.5),log10(1));
arSetPars('r8', log10(fixpars(12,1)),2,1,log10(1/16),log10(1/3));

arSetPars('rho', log10(fixpars(13,1)),2,1,log10(0.01),log10(0.9));
arSetPars('thet', log10(fixpars(14,1)),2,1,log10(0.01),log10(0.7));

arSetPars('und', log10(uu(1)),2,1,log10(0.002),log10(0.99));


%%%%

ar.config.optimizers=1;
%ar.config.showFitting = true;
arFitLHS(1)

arPlot;
close all

delete *mexa64


N0=ar.model.condition.xFineSimu(1,1);

time=ar.model.data.tFine;

qua=ar.model.data.yFineSimu(:,1);
hos=ar.model.data.yFineSimu(:,2);
icu=ar.model.data.yFineSimu(:,3);
rec=ar.model.data.yFineSimu(:,4);
dead=ar.model.data.yFineSimu(:,5);

datatoplot=[time,qua,hos,icu,rec,dead];

toplot(1:length(datatoplot),1:6)=datatoplot;

closevalue=[shift,ones(1,length(seq))];

se=[7:6:(length(seq)+1)*6];
params=10.^ar.p;
label=ar.pLabel;

R0=R0calc_IXD(params,label,uu(1),1);

savpars(:,1)=params';

arfixed=ar;

%% Reset before looping

params=10.^arfixed.p;

frac=arfixed.model.condition.xFineSimu(1,1)/N0;

R0=R0calc_IXD(params,label,uu(1),frac);
savpars(:,1)=params';
savchi2(:,1)=arfixed.chi2;
savr0(:,1)=R0;

clear ar;
ar=arfixed;


%% Moving window

k=0;
for i=seq
    disp(i)
    k=k+1;
    
    params=10.^ar.p;
    label=ar.pLabel;
    
    
    [minValue,closestIndex] = min(abs(ar.model.condition.tFine-closevalue(k)));
    
    %take the simulation's population at that index and save it in file for the
    %next fit (Models/last_model_t_loop_2.def)
    
    
    initcond = ar.model.condition.xFineSimu(closestIndex,1:15);
    savinit(:,k)=[params(find(strcmp(label,'r1'))),initcond]';
    rl=sprintf('init_x1     "%s"',initcond(1));
    
    Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',...
        cityname),...
        sprintf('Models/last_model_t_loop_3_%s.def',...
        cityname),'init_x1',rl)
    
    for j=2:15
        il=sprintf('init_x%d',j);
        rl=sprintf('init_x%d     "%d"',j,initcond(j));
        
        Func_replace_string(sprintf('Models/last_model_t_loop_3_%s.def',...
            cityname),...
            sprintf('Models/last_model_t_loop_3_%s.def',cityname),...
            il,rl)
    end
    
    
    [num,txt,raw]=xlsread(strcat('/Users/Arnab/Desktop/public_code/Testing_Model/',datatoread));
    minshift=0;
    for ii=i:i+datapoints-1;
        minshift=minshift+1;
        raw(ii,1)=mat2cell(cell2mat(raw(ii,1))-cell2mat(raw(ii,1))+minshift,1);
    end
    
    waslength=size(raw,1);
    dlrows=[2:i-1,i+datapoints:waslength];
    
    raw(dlrows,1:22)={[]};
    raw2=raw(setdiff(1:waslength,dlrows),:);
    %save the best dataset to fit
    
    filename=sprintf('Data/dat_%s_%d.xls',cityname,i);
    filenamed2d=sprintf('dat_%s_%d',cityname,i);
    xlwrite(filename, raw2);
    
    %start the job ...
    arInit;
    arLoadModel(sprintf('last_model_t_loop_3_%s',cityname));
    arLoadData(filenamed2d);
    arCompileAll;
    
    for ij=1:length(params)
        arSetPars(label{ij}, log10(params(ij)),2,1,-5,10);
    end
    
    %% parameters are fixed from Reference simulation
    % 'und' is varying
    %arSetPars('alpha', log10(fixpars(1,k+1)),2,1,log10(0.01),log10(0.4));
    arSetPars('bet', log10(fixpars(2,k+1)),2,1,log10(0.05),log10(0.1));
    arSetPars('delta', log10(fixpars(3,k+1)),2,1,log10(0.01),log10(0.9));
    
    %arSetPars('ki', log10(fixpars(4,k+1)),2,1,log10(10^-8),log10(1));
    
    arSetPars('r1', log10(fixpars(5,k+1)),2,1,log10(0.001),5);
    arSetPars('r10', log10(fixpars(6,k+1)),2,1,log10(0.1),log10(0.9));
    %arSetPars('r3', log10(fixpars(7,k+1)),2,1,log10(1/4.2),log10(2/5.2));
    %arSetPars('r4', log10(fixpars(8,k+1)),2,1,log10(1/14),log10(1/7));
    %arSetPars('r5', log10(fixpars(9,k+1)),2,1,log10(1/16),log10(1/5));
    %arSetPars('r6', log10(fixpars(10,k+1)),2,1,log10(1/7),log10(0.9));
    %arSetPars('r7', log10(fixpars(11,k+1)),2,1,log10(1/3.5),log10(1));
    %arSetPars('r8', log10(fixpars(12,k+1)),2,1,log10(1/16),log10(1/3));
    
    arSetPars('rho', log10(fixpars(13,k+1)),2,1,log10(0.01),log10(0.9));
    arSetPars('thet', log10(fixpars(14,k+1)),2,1,log10(0.01),log10(0.7));
    
    arSetPars('und', log10(uu(k)),2,1,-5,5);
    
    
    arFitLHS(1)
    
    arPlot
    
    time=ar.model.data.tFine;
    
    qua=ar.model.data.yFineSimu(:,1);
    hos=ar.model.data.yFineSimu(:,2);
    icu=ar.model.data.yFineSimu(:,3);
    rec=ar.model.data.yFineSimu(:,4);
    dead=ar.model.data.yFineSimu(:,5);
    
    datatoplot=[time,qua,hos,icu,rec,dead];
    toplot(1:length(datatoplot),se(k):se(k)+5)=datatoplot;
    
    arloop=ar;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% End
    params=10.^arloop.p;
    
    frac=arloop.model.condition.xFineSimu(1,1)/(N0-arloop.model.condition.xFineSimu(1,12));
    
    R0=R0calc_IXD(params,label,uu(1),frac);
    savpars(:,k+1)=params';
    savchi2(:,k+1)=arfixed.chi2;
    savr0(:,k+1)=R0;
    
    
    ar=arloop;
    
end


delete *mexmaci64

%% main folder needs to be adjusted as in 'fixed_undetected.m'
cd '/Users/Arnab/Desktop/public_code/Testing_Model/'

csvwrite(sprintf('area_spec_%s_fixed.csv',cityname),toplot)
save([sprintf('area_specific_%s_fixed.mat',cityname)],'savr0','savpars','savinit','savchi2','toplot','savr0per')


%     folder=pwd;

%     rmdir(strcat(cityname),'s');


end
