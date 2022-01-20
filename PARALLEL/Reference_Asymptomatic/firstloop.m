function ar=firstloop(datatoread,deffile,loopfile,cityname,missing,undetected,home_dir)

shift=5.2;
datapoints=7;
parstosave=15;

%     if undetected < 0.5
%         undetected=0.5;
%     end
%
mkdir(strcat(cityname))
mkdir(strcat(strcat(cityname),'/Models'))
mkdir(strcat(strcat(cityname),'/Data'))
copyfile(strcat('Models/',strcat(deffile,'.def')), strcat(strcat(cityname),'/Models'))
copyfile(strcat('Models/',loopfile), strcat(strcat(cityname),'/Models'))
copyfile('xlwrite.m',strcat(cityname))
copyfile('Copy_of_func_replace_string.m',strcat(cityname))
copyfile('R0calc.m',strcat(cityname))
%take the raw dataset
[num,txt,raw]=xlsread(datatoread);
cd(strcat(cityname))

seq=[3:length(raw)-datapoints+1];
%     seq=[3:5]; % for testing
susci=zeros(size(raw,1),1);
expo=zeros(size(raw,1),1);

savinit=zeros(15,length(seq)+1);
savpars=zeros(parstosave,length(seq)+1);
savchi2=zeros(1,length(seq)+1);
%%%adding the frac in savr0
savr0=zeros(2,length(seq)+1);
toplot=zeros(1000,6*(length(seq)+1));

savr0per=zeros(100,length(seq)+1);
savchi2per=zeros(100,length(seq)+1);
savfrac=zeros(100,length(seq)+1);

first_wind=15;
interval=first_wind-missing;
for i=2:interval+1
    raw(i,1)=mat2cell(cell2mat(raw(i,1))+shift-1,1);
end
raw(interval+2:end,1:22)={[]};

%save the new dataset to fit
%     xlwrite(sprintf('Data/dat_italy_d2d_tanmay_first_%s.xls',cityname), raw);
xlwrite(sprintf('Data/dat_first_%s.xls',cityname), raw);

%% Load D2D framework

%arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])

datafile=sprintf('dat_first_%s',cityname);

arInit;
arLoadModel(deffile);
arLoadData(datafile);

arCompileAll;

arSetPars('r1', log10(0.58),1,1,log10(0.001),log10(5));
arSetPars('r3', log10(1/4.2),1,1,log10(1/4.2),log10(2/5.2));
arSetPars('r4', log10(1/10),1,1,log10(1/14),log10(1/7));
%arSetPars('r4', log10(1/10),1,1,log10(1/14),log10(1/4));
arSetPars('r5', log10(1/10),1,1,log10(1/16),log10(1/5));
arSetPars('r6', log10(1/5),1,1,log10(1/7),log10(0.9));
arSetPars('r7', log10(1/2.5),1,1,log10(1/3.5),log10(1));
arSetPars('r8', log10(1/8),1,1,log10(1/16),log10(1/3));
arSetPars('r10', log10(0.5),1,1,log10(0.1),log10(0.9));

arSetPars('alpha', log10(0.4),2,1,log10(0.01),log10(0.4));
%arSetPars('bet', log10(1/4),1,1,log10(0.05),log10(1));
% arSetPars('bet', log10(0.05),1,1,log10(0.05),log10(0.25));
arSetPars('bet', log10(0.05),2,1,log10(0.05),log10(0.05));
% arSetPars('delta', log10(0.5),1,1,log10(0.3),log10(0.9));
arSetPars('delta', log10(0.5),1,1,log10(0.1),log10(0.9));
arSetPars('rho', log10(1/5),1,1,log10(0.01),log10(0.9));
arSetPars('thet', log10(0.2),1,1,log10(0.01),log10(0.7));

arSetPars('ki', log10(1),2,1,log10(10^-8),log10(1));
arSetPars('und', log10(undetected),2,1,log10(0.002),log10(0.99));

ar.config.optimizers=1;
% arFitLHS(50)
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
savr0(1,1)=R0;
savr0(2,1)=frac;

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
    initcond = ar.model.condition.xFineSimu(closestIndex,1:14);
    savinit(:,k)=[params(find(strcmp(label,'r1'))),initcond]';
    rl=sprintf('init_x1     "%s"',initcond(1));
    
    
    Copy_of_func_replace_string(sprintf('Models/model_loop_2_%s.def',...
        cityname),...
        sprintf('Models/model_loop_3_%s.def',...
        cityname),'init_x1',rl)
    
    for j=2:14
        il=sprintf('init_x%d',j);
        rl=sprintf('init_x%d     "%d"',j,initcond(j));
        %disp(rl)
        Copy_of_func_replace_string(sprintf('Models/model_loop_3_%s.def',...
            cityname),...
            sprintf('Models/model_loop_3_%s.def',cityname),...
            il,rl)
    end
    
    
    %     [num,txt,raw]=xlsread(strcat('/home/msi18/Desktop/covid-19/COVID/undetected/',datatoread));
    [num,txt,raw]=xlsread(datatoread);
    minshift=0;
    for ii=i:i+datapoints-1;
        minshift=minshift+1;
        raw(ii,1)=mat2cell(cell2mat(raw(ii,1))-cell2mat(raw(ii,1))+minshift,1);
    end
    
    dlrows=[2:i-1,i+datapoints:length(raw)];
    
    raw(dlrows,1:22)={[]};
    raw2=raw(setdiff(1:length(raw),dlrows),:);
    %save the best dataset to fit
    
    
    %         filename=sprintf('Data/dat_italy_d2d_tanmay_%s_%d.xls',cityname,i);
    %         filenamed2d=sprintf('dat_italy_d2d_tanmay_%s_%d',cityname,i);
    filename=sprintf('Data/dat_%s_%d.xls',cityname,i);
    filenamed2d=sprintf('dat_%s_%d',cityname,i);
    xlwrite(filename, raw2);
    
    
    %start the job ...
    arInit;
    arLoadModel(sprintf('model_loop_3_%s',cityname));
    arLoadData(filenamed2d);
    arCompileAll;
    
    for ij=1:length(params)
        arSetPars(label{ij}, log10(params(ij)),2,1,-5,10);
    end
    
    arSetPars('r1', log10(0.58),1,1,log10(0.001),log10(5));
    arSetPars('r10', log10(0.5),1,1,log10(0.1),log10(0.9));
    arSetPars('rho', log10(1/5),1,1,log10(0.01),log10(0.9));
    arSetPars('thet', log10(0.2),1,1,log10(0.01),log10(0.7));
    %     arSetPars('delta', log10(0.3),1,1,log10(0.3),log10(0.9));
    arSetPars('delta', log10(0.5),1,1,log10(0.1),log10(0.9));
    
    arFitLHS(20)
    
    arPlot
    
    time=ar.model.data.tFine;
    
    qua=ar.model.data.yFineSimu(:,1);
    hos=ar.model.data.yFineSimu(:,2);
    icu=ar.model.data.yFineSimu(:,3);
    rec=ar.model.data.yFineSimu(:,4);
    dead=ar.model.data.yFineSimu(:,5);
    
    datatoplot=[time,qua,hos,icu,rec,dead];
    toplot(1:length(datatoplot),se(k):se(k)+5)=datatoplot;
    
    %     close all
    
    %
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%% pertubation analysis Saham r0 sensitive vary
    %
    arloop=ar;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% End
    
    
    params=10.^arloop.p;
    
    frac=arloop.model.condition.xFineSimu(1,1)/(N0-arloop.model.condition.xFineSimu(1,12));
    R0=R0calc(params,label,frac);
    savpars(:,k+1)=params';
    savchi2(:,k+1)=arfixed.chi2;
    savr0(1,k+1)=R0;
    savr0(2,k+1)=frac;
    susci(k,1)=arloop.model.condition.xFineSimu(1,1);
    expo(k,1)=arloop.model.condition.xFineSimu(1,2);
    
    
    %         clear ar
    ar=arloop;
    
end
plabel=ar.pLabel;

delete Data/dat_italy_d2d_tanmay_*.*
delete *mexmaci64


cd(home_dir)% '/home/abp19/ITALYCOVID15122021/Ref_Asymp/'
rmdir(strcat(cityname),'s')

csvwrite(sprintf('area_spec_%s.csv',cityname),toplot)
save([sprintf('area_specific_%s.mat',cityname)],'savr0','savpars','savinit','savchi2','toplot','savr0per','plabel')
% csvwrite(sprintf('susc_%s.csv',cityname),susci)
% csvwrite(sprintf('expo_%s.csv',cityname),expo)

%     folder=pwd;

%     rmdir(strcat(cityname),'s');


end
