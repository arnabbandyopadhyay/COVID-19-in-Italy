function ar=fixedloop(datatoread,deffile,loopfile,cityname,missing,fixpars,pper)

shift=5.2; %incubation time
datapoints=7; %window size
parstosave=17; %parameters to save

%% create a folder with region name, copy all necessary files and change directory.
% all computation will happen in this folder
mkdir(strcat(cityname))

mkdir(strcat(strcat(cityname),'/Models'))
mkdir(strcat(strcat(cityname),'/Data'))
copyfile(strcat('Models/',strcat(deffile,'.def')), strcat(strcat(cityname),'/Models'))
copyfile(strcat('Models/',loopfile), strcat(strcat(cityname),'/Models'))
copyfile('xlwrite.m',strcat(cityname))
copyfile('Func_replace_string.m',strcat(cityname))
copyfile('R0calc.m',strcat(cityname))
%take the raw dataset
[num,txt,raw]=xlsread(datatoread);

%% take ICU and Hospital capacities
copyfile(strcat('icu_limits/',strcat(cityname,'_icu.xls')),strcat(cityname));
copyfile(strcat('hos_limits/',strcat(cityname,'_hos.xls')),strcat(cityname));

cd(strcat(cityname))

num2=xlsread(strcat(cityname,'_icu.xls'));
num3=xlsread(strcat(cityname,'_hos.xls'));

if missing>0
    vectlim=num2(missing:length(num2),2);
    hosplim=num3(missing:length(num3),2);
else
    vectlim=num2(1:length(num2),2);
    hosplim=num3(1:length(num3),2);
end

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

first_window=15;
interval=first_window;
for i=2:interval+1
    raw(i,1)=mat2cell(cell2mat(raw(i,1))+shift-1,1);
end
raw(interval+1:end,1:22)={[]};

%save the new dataset to fit
xlwrite(sprintf('Data/dat_first_%s.xls',cityname), raw);

%% Load D2D framework

%arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])

datafile=sprintf('dat_first_%s',cityname);

arInit;
arLoadModel(deffile);
arLoadData(datafile);

arCompileAll;

%% parameters are fixed as resulting from the Capacity model fit
% except for 'hlim' and 'icum': maximum capacities of hospital and icu, respectively
arSetPars('alpha', log10(fixpars(1,1)),2,1,-5,5);
arSetPars('bet', log10(fixpars(2,1)),2,1,-5,5);
arSetPars('delta', log10(fixpars(3,1)),2,1,-5,5);
arSetPars('hlim', log10(max(hosplim)),2,1,-10,10);
arSetPars('icum', log10(max(vectlim)),2,1,-10,10);
arSetPars('ki', log10(fixpars(6,1)),2,1,-5,5);
arSetPars('r1', log10(fixpars(7,1)),2,1,-5,5);
arSetPars('r10', log10(fixpars(8,1)),2,1,-5,5);
arSetPars('r3', log10(fixpars(9,1)),2,1,-5,5);
arSetPars('r4', log10(fixpars(10,1)),2,1,-5,5);
arSetPars('r5', log10(fixpars(11,1)),2,1,-5,5);
arSetPars('r6', log10(fixpars(12,1)),2,1,-5,5);
arSetPars('r7', log10(fixpars(13,1)),2,1,-5,5);
arSetPars('r8', log10(fixpars(14,1)),2,1,-5,5);
arSetPars('rho', log10(fixpars(15,1)),2,1,-5,5);
arSetPars('thet', log10(fixpars(16,1)),2,1,-5,5);
arSetPars('und', log10(fixpars(17,1)),2,1,-5,5);

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
R0=R0calc(params,label,1);

savpars(:,1)=params';

arfixed=ar;

%% Reset before looping

params=10.^arfixed.p;

frac=arfixed.model.condition.xFineSimu(1,1)/N0;
R0=R0calc(params,label,frac);
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
        %disp(rl)
        Func_replace_string(sprintf('Models/last_model_t_loop_3_%s.def',...
            cityname),...
            sprintf('Models/last_model_t_loop_3_%s.def',cityname),...
            il,rl)
    end
    
    
    [num,txt,raw]=xlsread(strcat('/Users/Arnab/Desktop/code/Cap_Model/',datatoread));
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
    
    arSetPars('alpha', log10(fixpars(1,k+1)),2,1,-5,5);
    arSetPars('bet', log10(fixpars(2,k+1)),2,1,-5,5);
    arSetPars('delta', log10(fixpars(3,k+1)),2,1,-5,5);
    arSetPars('hlim', log10(max(hosplim)),2,1,-10,10);
    arSetPars('icum', log10(max(vectlim)),2,1,-10,10);
    arSetPars('ki', log10(fixpars(6,k+1)),2,1,-5,5);
    arSetPars('r1', log10(fixpars(7,k+1)),2,1,-5,5);
    arSetPars('r10', log10(fixpars(8,k+1)),2,1,-5,5);
    arSetPars('r3', log10(fixpars(9,k+1)),2,1,-5,5);
    arSetPars('r4', log10(fixpars(10,k+1)),2,1,-5,5);
    arSetPars('r5', log10(fixpars(11,k+1)),2,1,-5,5);
    arSetPars('r6', log10(fixpars(12,k+1)),2,1,-5,5);
    arSetPars('r7', log10(fixpars(13,k+1)),2,1,-5,5);
    arSetPars('r8', log10(fixpars(14,k+1)),2,1,-5,5);
    arSetPars('rho', log10(fixpars(15,k+1)),2,1,-5,5);
    arSetPars('thet', log10(fixpars(16,k+1)),2,1,-5,5);
    arSetPars('und', log10(fixpars(17,k+1)),2,1,-5,5);
    
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
    R0=R0calc(params,label,frac);
    savpars(:,k+1)=params';
    savchi2(:,k+1)=arfixed.chi2;
    savr0(:,k+1)=R0;
    
    ar=arloop;
    
end

delete *mexmaci64


%% main folder needs to be adjusted as in 'main_fixed.m'
cd '/Users/Arnab/Desktop/code/Cap_Model/'

csvwrite(sprintf('capacity_%s_fixed_%d.csv',cityname,pper),toplot)
save([sprintf('capacity_%s_fixed_%d.mat',cityname,pper)],'savr0','savpars','savinit','savchi2','toplot','savr0per')


%     folder=pwd;

%     rmdir(strcat(cityname),'s');


end
