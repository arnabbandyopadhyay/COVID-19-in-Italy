function ar=perturb_loop(datatoread,deffile,loopfile,cityname,missing,undetected)%,para)

    shift=5.2; % incubation period
    datapoints=7; % window size
    npertu=1000; % perturbation trial
    parstosave=15; % parameters to save
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
    cd(strcat(cityname))
    
    %% Arrays to store relevant information
    
    seq=[3:length(raw)-datapoints+1];
%     seq=[3:5]; % for testing

    savinit=zeros(15,length(seq)+1);
    savpars=zeros(parstosave,length(seq)+1);
    savchi2=zeros(1,length(seq)+1);
    savr0=zeros(1,length(seq)+1);
    toplot=zeros(1000,6*(length(seq)+1));

    savr0per=zeros(npertu,length(seq));
    savchi2per=zeros(100,length(seq)+1);
    savfrac=zeros(100,length(seq)+1);

    paramsperturb=zeros(npertu,length(seq)+1);
    
    %% Creating data frame by taking first 15 days = up to 9th March
    
    first_wind=15;
    interval=(first_wind-missing);
    
    for i=2:interval+1
        raw(i,1)=mat2cell(cell2mat(raw(i,1))+shift-1,1);
    end
    raw(interval+1:end,1:22)={[]};
    
    %save the new dataset to fit
    xlwrite(sprintf('Data/dat_italy_d2d_first_%s.xls',cityname), raw);

    %% Load D2D framework  
    
    datafile=sprintf('dat_italy_d2d_first_%s',cityname);
    
    arInit; % Initializing D2D framework
    arLoadModel(deffile); % loading model
    arLoadData(datafile); % loading data

    arCompileAll; % compiling data with model
    
    % Declaration of parameters, their starting value, searching range
    %arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])  
    
    arSetPars('r1', log10(0.58),1,1,log10(0.001),5);

    arSetPars('r4', log10(1/10),1,1,log10(1/14),log10(1/7));
    
    arSetPars('r5', log10(1/10),1,1,log10(1/16),log10(1/5));
    arSetPars('r6', log10(1/5),1,1,log10(1/7),log10(0.9));
    
    arSetPars('r7', log10(1/2.5),1,1,log10(1/3.5),log10(1));
    arSetPars('r8', log10(1/8),1,1,log10(1/16),log10(1/3));
    arSetPars('r10', log10(0.5),1,1,log10(0.3),log10(0.9));
    
    arSetPars('rho', log10(1/5),1,1,log10(0.01),log10(0.9));
    
    arSetPars('thet', log10(0.2),1,1,log10(0.01),log10(0.7));
    
    arSetPars('alpha', log10(0.4),2,1,log10(0.01),log10(0.4));
    
    arSetPars('bet', log10(1/4),1,1,log10(0.05),log10(1));
    
    arSetPars('delta', log10(0.15),1,1,log10(0.0001),log10(0.9));
    arSetPars('r3', log10(1/4.2),1,1,log10(1/4.2),log10(2/5.2));
    
    arSetPars('ki', log10(1),2,1,log10(10^-8),log10(1));
    arSetPars('und', log10(undetected),2,1,log10(0.002),log10(0.99));

    ar.config.optimizers=1;
    %ar.config.showFitting = true;
%     arFitLHS(50)
    arFitLHS(10)
    arPlot;
    close all
    
    delete *mexa64
    
    %% storing information to plot after fit 
    
    N0=ar.model.condition.xFineSimu(1,1); % Initial population

    time=ar.model.data.tFine;

    qua=ar.model.data.yFineSimu(:,1);   % infection
    hos=ar.model.data.yFineSimu(:,2);   % hospital
    icu=ar.model.data.yFineSimu(:,3);   % ICU
    rec=ar.model.data.yFineSimu(:,4);   % recovery
    dead=ar.model.data.yFineSimu(:,5);  % dead

    
    datatoplot=[time,qua,hos,icu,rec,dead];

    toplot(1:length(datatoplot),1:6)=datatoplot;

    closevalue=[shift,ones(1,length(seq))];

    se=[7:6:(length(seq)+1)*6];
    
    %% calculating Rt
    
    params=10.^ar.p; % parameters are stored in ar.p variable in log scale
    label=ar.pLabel; % corresponding name of parameter
    R0=R0calc(params,label,1);  % calculating Rt
    
    savpars(:,1)=params';   % 
    
    
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
    
    plbper=params-params.*0.2;
    pubper=params+params.*0.2;
    
    % upper and lower bound of parameters
    
    alphar=[0.01,0.4];
    betar=[0.05,1];
    r10r=[0.3,0.9];
    deltar=[0.1,0.9];
    kir=[1,1];
    r1r=[0,1];
    r3r=[1/4.2,2/5.2];
    r4r=[1/14,1/7];
    r5r=[1/16,1/5];
    r6r=[1/7,0.9];
    r7r=[1/3.5,1];
    r8r=[1/16,1/3];
    rhor=[0.1,1];
    thetr=[0.1,1];
    

    lbper=[alphar(1),betar(1),r10r(1),deltar(1),kir(1),r1r(1),r3r(1),r4r(1),r5r(1),r6r(1),r7r(1),r8r(1),rhor(1),thetr(1)];
    ubper=[alphar(2),betar(2),r10r(2),deltar(2),kir(2),r1r(2),r3r(2),r4r(2),r5r(2),r6r(2),r7r(2),r8r(2),rhor(2),thetr(2)];
    
    % if perturbing parameter exceeds upper and lower bound of the parameter,
    % reset it to its own upper and lower bound
    
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
    for perturb=1:perturb_todo % length(lty)
        
        ar=arfixed;
        k=0;
        for i=seq
            disp(i)
            fprintf('no of per=%d\n',length(lty));
            fprintf('loop per=%d\n',perturb);
            k=k+1;
    %         undetected=uu(k);

%             params=10.^ar.p;
            label=ar.pLabel;


            [minValue,closestIndex] = min(abs(ar.model.condition.tFine-closevalue(k)));

            % store the state variables for the next window
            
            initcond = ar.model.condition.xFineSimu(closestIndex,1:14);
            savinit(:,k)=[params(find(strcmp(label,'r1'))),initcond]';
            rl=sprintf('init_x1     "%d"',initcond(1));


            Func_replace_string(sprintf('Models/last_model_t_loop_2_%s.def',...
                cityname),...
                sprintf('Models/last_model_t_loop_3_%s.def',...
                cityname),'init_x1',rl)

            for j=2:14
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

            [num,txt,raw]=xlsread(strcat('/Users/Arnab/Desktop/code/Ref_and_Asym_Model/',datatoread));
            minshift=0;
            for ii=i:i+datapoints-1;
                minshift=minshift+1;
                raw(ii,1)=mat2cell(cell2mat(raw(ii,1))-cell2mat(raw(ii,1))+minshift,1);
            end

            dlrows=[2:i-1,i+datapoints:length(raw)];

            raw(dlrows,1:22)={[]};
            raw2=raw(setdiff(1:length(raw),dlrows),:);
            %save the best dataset to fit

            % load model, data and start window fitting
            
            filename=sprintf('Data/dat_italy_d2d_%s_%d.xls',cityname,i);
            filenamed2d=sprintf('dat_italy_d2d_%s_%d',cityname,i);
            xlwrite(filename, raw2);

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
            
            % fractional parameters and R1, R10 will adapt in each window
            
            arSetPars('r1', log10(0.58),1,1,log10(0.001),log10(3));
            arSetPars('rho', log10(1/5),1,1,log10(0.01),log10(0.9));
            arSetPars('thet', log10(0.2),1,1,log10(0.01),log10(0.7));
            arSetPars('delta', log10(0.15),1,1,log10(0.01),log10(0.9));
            arSetPars('r10', log10(0.5),1,1,log10(0.3),log10(0.9));
            arSetPars('und', log10(undetected),2,1,log10(0.002),log10(0.99));
           
            arFitLHS(10)

            arPlot
            close all

            % save observations in each window
            
            time=ar.model.data.tFine;

            qua=ar.model.data.yFineSimu(:,1); % infection
            hos=ar.model.data.yFineSimu(:,2);   % hospital
            icu=ar.model.data.yFineSimu(:,3);   % ICU
            rec=ar.model.data.yFineSimu(:,4);   % recovery
            dead=ar.model.data.yFineSimu(:,5);  % dead

            datatoplot=[time,qua,hos,icu,rec,dead];
            toplot(1:length(datatoplot),se(k):se(k)+5)=datatoplot;

        %     close all

            arloop=ar;

            % calculate Rt in each window 
            
            params=10.^arloop.p;
            
            frac=arloop.model.condition.xFineSimu(1,1)/(N0-arloop.model.condition.xFineSimu(1,12));
            R0=R0calc(params,label,frac);
            disp([params(1:end),label,N0,R0])
            savpars(:,k+1)=params';
            savchi2(:,k+1)=arfixed.chi2;
            savr0(:,k+1)=R0;
%             
            savr0per(perturb,k)=R0;
%             paramsperturb(perturb,k)=params;

    %         clear ar
            ar=arloop;

        end
        TP(perturb).toplot=toplot;
        TP(perturb).params=savpars;
        TP(perturb).r0=savr0per;
        
        delete Data/dat_italy_d2d_*.*
        delete *mexmaci64
        
    end
    
    
    %% main path needs to be adjusted as in final_model_perturb.m, change accordingly
    cd '/Users/Arnab/Desktop/code/Ref_and_Asym_Model/'
    
%     csvwrite(sprintf('area_spec_%s.csv',cityname),toplot)
    save([sprintf('perturbation_%s.mat',cityname)],'savr0per','TP')
%     
    
%     folder=pwd;
    
%     rmdir(strcat(cityname),'s');
    csvwrite(sprintf('final_r0_%s.csv',cityname),TP(perturb_todo).r0(1:perturb_todo,:));
    
    
    
end


    
