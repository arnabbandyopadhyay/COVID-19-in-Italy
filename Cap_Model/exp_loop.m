function ar=exp_loop(datatoread,deffile,loopfile,cityname,missing,undetected)

shift=5.2; % incubation period
npertu=200;
datapoints=7; % window size
parstosave=17; % parameters to save
perturb_todo=30; % perturbation

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

copyfile(strcat('icu_limits/',strcat(cityname,'_icu.xls')),strcat(cityname));
copyfile(strcat('hos_limits/',strcat(cityname,'_hos.xls')),strcat(cityname));

cd(strcat(cityname))

%% Create vectors for hospital and ICU capacities in the windows
num2=xlsread(strcat(cityname,'_icu.xls'));
num3=xlsread(strcat(cityname,'_hos.xls'));

if missing>0
    vectlim=num2((missing+1):length(num2),2);
    hosplim=num3((missing+1):length(num3),2);
else
    vectlim=num2(1:length(num2),2);
    hosplim=num3(1:length(num3),2);
end

seq=[3:length(raw)-datapoints+1];
%     seq=[3:5]; % for testing

%% Arrays to store relevant information
savinit=zeros(16,length(seq)+1);

savpars=zeros(parstosave,length(seq)+1);
savchi2=zeros(1,length(seq)+1);
savr0=zeros(1,length(seq)+1);
toplot=zeros(1000,6*(length(seq)+1));

savr0per=zeros(npertu,length(seq)+1);
savchi2per=zeros(100,length(seq)+1);
savfrac=zeros(100,length(seq)+1);

paramsperturb=zeros(npertu,length(seq)+1);

%% Creating data frame by taking the first 15 days = up to 9th March
if missing>0
    first_window=15-missing;
else
    first_window=15;
end


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

arSetPars('r1', log10(0.58),1,1,log10(0.001),5);
arSetPars('r3', log10(1/4.2),1,1,log10(1/4.2),log10(2/5.2));
arSetPars('r4', log10(1/10),1,1,log10(1/14),log10(1/7));
%     arSetPars('r4', log10(1/10),1,1,log10(1/14),log10(1/4));
arSetPars('r5', log10(1/10),1,1,log10(1/16),log10(1/5));
arSetPars('r6', log10(1/5),1,1,log10(1/7),log10(0.9));
arSetPars('r7', log10(1/2.5),1,1,log10(1/3.5),log10(1));
arSetPars('r8', log10(1/8),1,1,log10(1/16),log10(1/3));
arSetPars('r10', log10(0.5),1,1,log10(0.1),log10(0.9));

arSetPars('alpha', log10(0.4),2,1,log10(0.01),log10(0.4));
arSetPars('bet', log10(1/4),1,1,log10(0.05),log10(1));
arSetPars('delta', log10(0.5),1,1,log10(0.3),log10(0.9));
arSetPars('rho', log10(1/5),1,1,log10(0.01),log10(0.9));
arSetPars('thet', log10(0.2),1,1,log10(0.01),log10(0.7));

arSetPars('ki', log10(1),2,1,log10(10^-8),log10(1));
arSetPars('und', log10(undetected),2,1,log10(0.002),log10(0.99));

%     arSetPars('hlim', log10(max(hosplim)),2,1,-10,10);
%     arSetPars('icum', log10(max(vectlim)),2,1,-10,10);

arSetPars('icum', log10(vectlim(1)),2,1,-5,5);
arSetPars('hlim', log10(hosplim(1)),2,1,-5,5);

ar.config.optimizers=1;

arFitLHS(50)


arPlot;
close all

delete *mexa64

%% storing information to plot after fit

N0=ar.model.condition.xFineSimu(1,1); % Initial population

time=ar.model.data.tFine;

qua=ar.model.data.yFineSimu(:,1); % infection
hos=ar.model.data.yFineSimu(:,2); % hospital
icu=ar.model.data.yFineSimu(:,3); % ICU
rec=ar.model.data.yFineSimu(:,4); % recovered
dead=ar.model.data.yFineSimu(:,5); % dead


datatoplot=[time,qua,hos,icu,rec,dead];

toplot(1:length(datatoplot),1:6)=datatoplot;

closevalue=[shift,ones(1,length(seq))];

se=[7:6:(length(seq)+1)*6];

%% calculating Rt
params=10.^ar.p; % parameters are stored in ar.p variable in log scale
label=ar.pLabel; % corresponding name of parameters
R0=R0calc(params,label,1); % calculating Rt

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

%% perturbation 20% of its base value

caplb=params-params.*0.2;
capub=params+params.*0.2;

pid=[1:3,6:16];

plbper=caplb(pid);
pubper=capub(pid);

%% upper and lower bound of parameters

alphar=[0.01,0.4];
betar=[0.05,1];
deltar=[0.1,0.9];
kir=[1,1];
r1r=[0,1];
r10r=[0.3,0.9];
r3r=[1/4.2,2/5.2];
r4r=[1/14,1/7];
r5r=[1/16,1/5];
r6r=[1/7,0.9];
r7r=[1/3.5,1];
r8r=[1/16,1/3];
rhor=[0.1,1];
thetr=[0.1,1];


lbper=[alphar(1),betar(1),deltar(1),kir(1),r1r(1),r10r(1),r3r(1),r4r(1),r5r(1),r6r(1),r7r(1),r8r(1),rhor(1),thetr(1)];
ubper=[alphar(2),betar(2),deltar(2),kir(2),r1r(2),r10r(2),r3r(2),r4r(2),r5r(2),r6r(2),r7r(2),r8r(2),rhor(2),thetr(2)];

for nn=1:length(lbper)
    
    if plbper(nn)<lbper(nn)
        plbper(nn)=lbper(nn);
    end
    
    if pubper(nn)>ubper(nn)
        pubper(nn)=ubper(nn);
    end
    
end

% once upper and lower bound set properly, sample npertu times
% within the range and find which combination of parameters
% results 20% or less variation in terms of total parameter variation
% (Barkai and Leibler, Nature, 1997). once we have perturbed parameter
% set, we proceed to moving window.

alphaper=random('uniform',plbper(1),pubper(1),[npertu,1]);
betper=random('uniform',plbper(2),pubper(2),[npertu,1]);
kiper=random('uniform',1,1,[npertu,1]);
r3per=random('uniform',plbper(7),pubper(7),[npertu,1]);
r4per=random('uniform',plbper(8),pubper(8),[npertu,1]);
r5per=random('uniform',plbper(9),pubper(9),[npertu,1]);
r6per=random('uniform',plbper(10),pubper(10),[npertu,1]);
r7per=random('uniform',plbper(11),pubper(11),[npertu,1]);
r8per=random('uniform',plbper(12),pubper(12),[npertu,1]);

alpha0=params(find(strcmp(label,'alpha')));
bet0=params(find(strcmp(label,'bet')));
ki0=params(find(strcmp(label,'ki')));
r30=params(find(strcmp(label,'r3')));
r40=params(find(strcmp(label,'r4')));
r50=params(find(strcmp(label,'r5')));
r60=params(find(strcmp(label,'r6')));
r70=params(find(strcmp(label,'r7')));
r80=params(find(strcmp(label,'r8')));

k0=alpha0*bet0*ki0*r30*r40*r50*r60*r70*r80;
lty=[];
for tt=1:npertu
    kperturb=alphaper(tt)*betper(tt)*kiper(tt)*r3per(tt)*r4per(tt)*r5per(tt)*r6per(tt)*r7per(tt)*r8per(tt);
    if abs(1-kperturb/k0) <= 0.2
        lty=[lty,tt];
        %         else
        %             lty=[lty,kperturb/k0];
        %
    end
    
end




%% Moving window

for perturb=1:perturb_todo %length(lty)
    ar=arfixed;
    
    k=0;
    toplot(:,7:(6*(length(seq)+1)))=0;
    savpars(:,2:(length(seq)+1))=0;
    for i=seq
        disp(i)
        fprintf('no of per=%d\n',length(lty));
        fprintf('loop per=%d\n',perturb);
        k=k+1;
        
        params=10.^ar.p;
        label=ar.pLabel;
        
        
        [minValue,closestIndex] = min(abs(ar.model.condition.tFine-closevalue(k)));
        
        % store the state variables for the next window
        
        initcond = ar.model.condition.xFineSimu(closestIndex,1:15);
        savinit(:,k)=[params(find(strcmp(label,'r1'))),initcond]';
        rl=sprintf('init_x1     "%s"',initcond(1));
        
        modelfile=sprintf('Models/last_model_t_loop_3_%s.def',cityname);
        
        if isfile(modelfile)
            delete(modelfile);
            Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',...
                cityname),...
                sprintf('Models/last_model_t_loop_3_%s.def',...
                cityname),'init_x1',rl)
            
        else
            Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',...
                cityname),...
                sprintf('Models/last_model_t_loop_3_%s.def',...
                cityname),'init_x1',rl)
        end
        
        for j=2:15
            il=sprintf('init_x%d',j);
            rl=sprintf('init_x%d     "%d"',j,initcond(j));
            %disp(rl)
            Func_replace_string(sprintf('Models/last_model_t_loop_3_%s.def',...
                cityname),...
                sprintf('Models/last_model_t_loop_3_%s.def',cityname),...
                il,rl)
        end
        
        % load the original data and frame shift accordingly, change
        % the path accordingly
        
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
        % load model, data and start window fitting
        
        filename=sprintf('Data/dat_%s_%d.xls',cityname,i);
        filenamed2d=sprintf('dat_%s_%d',cityname,i);
        
        if isfile(filename)
            delete(filename);
            xlwrite(filename, raw2);
        else
            xlwrite(filename, raw2);
        end
        
        %start the job ...
        arInit;
        arLoadModel(sprintf('last_model_t_loop_3_%s',cityname));
        arLoadData(filenamed2d);
        arCompileAll;
        
        % declare parameters values, physiological parameters will be
        % fixed throughout
        
        arSetPars('alpha', log10(alphaper(lty(perturb))),2,1,-5,5);
        arSetPars('bet', log10(betper(lty(perturb))),2,1,-5,5);
        arSetPars('ki', log10(kiper(lty(perturb))),2,1,-5,5);
        arSetPars('r3', log10(r3per(lty(perturb))),2,1,-5,5);
        arSetPars('r4', log10(r4per(lty(perturb))),2,1,-5,5);
        arSetPars('r5', log10(r5per(lty(perturb))),2,1,-5,5);
        arSetPars('r6', log10(r6per(lty(perturb))),2,1,-5,5);
        arSetPars('r7', log10(r7per(lty(perturb))),2,1,-5,5);
        arSetPars('r8', log10(r8per(lty(perturb))),2,1,-5,5);
        
        
        arSetPars('r1', log10(0.58),1,1,log10(0.001),5);
        arSetPars('r10', log10(0.5),1,1,log10(0.1),log10(0.9));
        arSetPars('rho', log10(1/5),1,1,log10(0.01),log10(0.9));
        arSetPars('thet', log10(0.2),1,1,log10(0.01),log10(0.7));
        arSetPars('delta', log10(0.15),1,1,log10(0.01),log10(0.9));
        
        arSetPars('icum', log10(vectlim(k)),2,1,-5,5);
        arSetPars('hlim', log10(hosplim(k)),2,1,-5,5);
        arSetPars('und', log10(undetected),2,1,log10(0.002),log10(0.99));
        %%
        
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
        
        close all
        arloop=ar;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% End
        
        %% calculate Rt in each window
        
        params=10.^arloop.p;
        
        frac=arloop.model.condition.xFineSimu(1,1)/(N0-arloop.model.condition.xFineSimu(1,12));
        R0=R0calc(params,label,frac);
        savpars(:,k+1)=params';
        savchi2(:,k+1)=arfixed.chi2;
        savr0(:,k+1)=R0;
        
        %         clear ar
        ar=arloop;
        
    end
    
    TP(perturb).toplot=toplot;
    TP(perturb).params=savpars;
    TP(perturb).r0=savr0per;
    
end

delete *mexmaci64


%% main path needs to be adjusted as in final_model_perturb.m, change accordingly
cd '/Users/Arnab/Desktop/code/Cap_Model/'

save([sprintf('area_specific_%s.mat',cityname)],'savr0per','TP')

csvwrite(sprintf('area_spec_%s.csv',cityname),toplot)
save([sprintf('parfit_%s.mat',cityname)],'savr0','savpars','savinit','savchi2','toplot','savr0per')


end
