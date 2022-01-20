% close all
clear all
clc

res_toPerturb="CapacResults";

num_pert=100;
var_toApply=0.1;

mkdir('cap_res')
mkdir('deCap_res')
mkdir('EarlyTest_res')
mkdir('LateTest_res')

mkdir('dead_diff')

mkdir('rt_cap')


mainfol=pwd;

regions = [{'lombardia'}, {'abruzzo'}, {'campania'}, {'emilia'}, {'friuli'},...
    {'lazio'}, {'liguria'}, {'marche'}, ...
    {'molise'}, {'piemonte'}, {'sardegna'}, {'sicilia'}, {'toscana'},...
    {'valledaosta'}, {'veneto'}, {'calabria'}, {'puglia'}, {'umbria'},...
    {'basilicata'}, {'bolzano'}, {'trento'}];
% regions = [{'lombardia'}];
%{'italy'},
% regions=[{'lombardia'}, {'emilia'}, {'liguria'}, {'marche'}, ...
%     {'piemonte'}, {'valledaosta'}, {'bolzano'}];

P0.abruzzo=1311580;P0.campania=5801692;P0.emilia=4459477;P0.friuli=1215220;
P0.lazio=5879082;P0.liguria=1550640;P0.lombardia=10060574;P0.italy=60359546;
P0.marche=1525271;P0.molise=305617;P0.piemonte=4356406;P0.sardegna=1639591;
P0.sicilia=4999891;P0.toscana=3729641;P0.valledaosta=125666;P0.veneto=4905854;
P0.calabria=1947131;P0.puglia=4029053;P0.umbria=882015;P0.basilicata=562869;
P0.bolzano=520891;P0.trento=538223;


%%%initial conditions: undetected model
C0.abruzzo=1*1/(1-0.8939); C0.campania=3*1/(1-0.1856); C0.emilia=18*1/(1-0.9169); C0.friuli=6*1/(1-0.891);
C0.italy=221*1/(1-0.9357); C0.lazio=3*1/(1-0.5); C0.liguria=1*1/(1-0.9309); C0.lombardia=225*1/(1-0.9553);
C0.marche=1*1/(1-0.9057); C0.molise=3*1/(1-0.5186); C0.piemonte=2*1/(1-0.9294); C0.sardegna=1*1/(1-0.9553);
C0.sicilia=1*1/(1-0.5); C0.toscana=2*1/(1-0.8777); C0.valledaosta=2*1/(1-0.9025); C0.veneto=32*1/(1-0.8967);
C0.calabria=1*1/(1-0.8442); C0.puglia=1*1/(1-0.9433); C0.umbria=2*1/(1-0.5); C0.basilicata=1*1/(1-0.5733);
C0.bolzano=1*1/(1-0.9261); C0.trento=4*1/(1-0.9161);

CN.abruzzo='Abruzzo'; CN.campania='Campania'; CN.emilia='Emilia'; CN.friuli='Friuli';
CN.italy='Italy'; CN.lazio='Lazio'; CN.liguria='Liguria'; CN.lombardia='Lombardia';
CN.marche='Marche'; CN.molise='Molise'; CN.piemonte='Piemonte'; CN.sardegna='Sardegna';
CN.sicilia='Sicilia'; CN.toscana='Toscana'; CN.valledaosta='Valledaosta'; CN.veneto='Veneto';
CN.calabria='Calabria'; CN.puglia='Puglia'; CN.umbria='Umbria'; CN.basilicata='Basilicata';
CN.bolzano='Bolzano'; CN.trento='Trento';

%%%%%%%
mm=1;
getthefile1=strcat(res_toPerturb,'/area_specific_',getfield(CN,regions{mm}),'.mat');
fileit=load(getthefile1);

disp(getfield(CN,regions{mm}))
pop=getfield(P0,regions{mm});
param=fileit.savpars;
maxtime=(size(param,2));
if mm==1 refval=maxtime; missval=0; else missval=refval-maxtime; end
%%%%%%

parfor mm=1:length(regions)
    DECAP=1;
    
    %%
    %load
    getthefile1=strcat(res_toPerturb,'/area_specific_',getfield(CN,regions{mm}),'.mat');
    fileit=load(getthefile1);
    
    disp(getfield(CN,regions{mm}))
    pop=getfield(P0,regions{mm});
    %disp(pop)
    %%
    %save parameters and initial conditions
    
    param=fileit.savpars;
    initco=fileit.savinit(2:size(fileit.savinit,1),:);
    %             maxtime=(size(param,2)-1);
    maxtime=(size(param,2));
    missval=refval-maxtime;
    disp(maxtime)
    %%
    %how to proceed:
    %A:save initial conditions and parameters from the first window.
    
    CHOOSEINI=1;
    
    ini_par=param(:,CHOOSEINI);
    ini_init=initco(:,1);
    
    r1r=[(0.001),5];%7
    r3r=[(1/4.2),(2/5.2)];%8
    r4r=[(1/14),(1/7)];%9
    
    r5r=[(1/16),(1/5)];%10
    r6r=[(1/7),(0.9)];%11
    r7r=[(1/3.5),(1)];%12
    r8r=[(1/16),(1/3)];%13
    r10r=[(0.1),(0.9)];%6
    
    alphar=[(0.01),(0.5)];%1
    
    betr=[(0.05),(0.1)];%2
    deltar=[(0.1),(0.9)];%3
    rhor=[(0.01),(0.9)];%14
    thetr=[(0.01),(0.7)];%15
    
    kir=[(10^-8),(1)];%6
    undr=[(0.002),(0.99)];%16
    
    alpha_0=ini_par(1);
    bet_0=ini_par(2);
    delta_0=ini_par(3);
    
    ki_0=ini_par(6);
    r1_0=ini_par(7);
    r10_0=ini_par(8);
    r3_0=ini_par(9);
    r4_0=ini_par(10);
    r5_0=ini_par(11);
    r6_0=ini_par(12);
    r7_0=ini_par(13);
    r8_0=ini_par(14);
    rho_0=ini_par(15);
    thet_0=ini_par(16);
    und=ini_par(17);
    
    rhos=param(15,CHOOSEINI:maxtime);
    thets=param(16,CHOOSEINI:maxtime);
    r1s=param(7,CHOOSEINI:maxtime);
    r10s=param(8,CHOOSEINI:maxtime);
    deltas=param(3,CHOOSEINI:maxtime);
    
    lbper=[alphar(1),betr(1),deltar(1),kir(1),r1r(1),r10r(1),r3r(1),...
        r4r(1),r5r(1),r6r(1),r7r(1),r8r(1),rhor(1),thetr(1)];
    ubper=[alphar(2),betr(2),deltar(2),kir(2),r1r(2),r10r(2),r3r(2),...
        r4r(2),r5r(2),r6r(2),r7r(2),r8r(2),rhor(2),thetr(2)];
    
    lbper2=[rhor(1),thetr(1),r10r(1),deltar(1),r1r(1)];
    ubper2=[rhor(2),thetr(2),r10r(2),deltar(2),r1r(2)];
    
    
    %%
    %B:apply +/-10% to the fractions (und excluded) WITHIN their range.
    rhoperF=[]; thetperF=[]; r10perF=[]; deltaperF=[]; r1perF=[];
    for outind=CHOOSEINI:maxtime
        
        k0=rhos(outind)*thets(outind)*r10s(outind)*deltas(outind)*r1s(outind);
        base_par=[rhos(outind),thets(outind),r10s(outind),deltas(outind),r1s(outind)];
        
        ini_par=param;%(:,2); ???MS
        
        plbper=base_par-base_par.*var_toApply;
        pubper=base_par+base_par.*var_toApply;
        for parlist=1:size(plbper,1)
            if plbper(parlist)<lbper2(parlist)
                plbper(parlist)=lbper2(parlist);
                disp('lower')
            end
            if pubper(parlist)>ubper2(parlist)
                pubper(parlist)=ubper2(parlist);
                disp('upper')
            end
        end
        %%
        rhoper=[]; thetper=[]; r10per=[]; deltaper=[]; r1per=[];
        %%
        %C:save parameters combinations.
        sizelt=0;
        
        while sizelt < num_pert
            
            rhoper1=random('uniform',plbper(1),pubper(1),[1,1]);
            thetper1=random('uniform',plbper(2),pubper(2),[1,1]);
            r10per1=random('uniform',plbper(3),pubper(3),[1,1]);
            deltaper1=random('uniform',plbper(4),pubper(4),[1,1]);
            r1per1=random('uniform',plbper(5),pubper(5),[1,1]);
            
            kperturb=rhoper1*thetper1*r10per1*deltaper1*r1per1;
            if abs(1-kperturb/k0) <= var_toApply
                
                rhoper=[rhoper,rhoper1];
                thetper=[thetper,thetper1];
                r10per=[r10per,r10per1];
                deltaper=[deltaper,deltaper1];
                r1per=[r1per,r1per1];
                
            end
            
            sizelt=length(rhoper);
            
        end
        
        rhoperF=[rhoperF,rhoper'];
        thetperF=[thetperF,thetper'];
        r10perF=[r10perF,r10per'];
        deltaperF=[deltaperF,deltaper'];
        r1perF=[r1perF,r1per'];
        
    end
    rhoperF=[rhoperF;rhos];
    thetperF=[thetperF;thets];
    r10perF=[r10perF;r10s];
    deltaperF=[deltaperF;deltas];
    r1perF=[r1perF;r1s];
    %D:run each; save; plot.
    
    num2=xlsread(strcat('../icu_limits/',getfield(CN,regions{mm}),'_icu.xls'));
    num3=xlsread(strcat('../hos_limits/',getfield(CN,regions{mm}),'_hos.xls'));
    disp('missing')
    disp(missval)
    %     if missval>0
    %         vectlim=num2((missval+1):length(num2),2);
    %         hosplim=num3((missval+1):length(num3),2);
    %     else
    %         vectlim=num2(1:length(num2),2);
    %         hosplim=num3(1:length(num3),2);
    %     end
    %%%%%%%
    hosplim=param(4,:);
    vectlim=param(5,:);
    %%%%%%%%%%%
    for p_num=1:(num_pert+1)
        disp('time')
        disp(maxtime)
        tsave=[]; xsave=[]; rtsave=[]; parsave=[]; x0save=[]; deadsave=[];
        for tt=1:(maxtime-1)
            if tt==1
                tspan=0:0.2:5.2;
            elseif tt<(maxtime-1)
                tspan=0:1;%16July
            else
                tspan=0:6;
            end
            rho_0=rhoperF(p_num,tt);
            thet_0=thetperF(p_num,tt);
            r10_0=r10perF(p_num,tt);
            delta_0=deltaperF(p_num,tt);
            r1_0=r1perF(p_num,tt);
            hlim=hosplim(tt);
            icum=vectlim(tt);
            if tt==1
                %                         x0=ini_init;
                x0=[getfield(P0,regions{mm}) getfield(C0,regions{mm}) repelem(0,(15-2))];
            else
                x0=xsave(tt-1,:);
            end
            %             dbstop if naninf
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
            [t,x]=ode15s(@(t,x) ODEcovCAP(t,x, alpha_0, bet_0, delta_0, ki_0,...
                r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                und, pop, hlim, icum),tspan,x0);%,opts);
            
            %%calc Rt
            r9=2*r3_0*r4_0/(2*r4_0+r3_0);
            m=((1-und)/(1-alpha_0));
            num=r1_0*(alpha_0/r9 + (1-alpha_0)/r3_0 + bet_0*(1-alpha_0)*(1-rho_0)*m/r4_0...
                +bet_0*rho_0*(1-alpha_0)*m/r6_0 + ki_0*(1-alpha_0)*(1-m)/r4_0);
            frac=x(size(x,1),1)/(pop-x(size(x,1),12));
            rtsave=[rtsave; num*frac];
            parsave=[parsave;[alpha_0 bet_0 delta_0 ki_0 r1_0 r10_0 r3_0 r4_0 r5_0 r6_0 r7_0 r8_0 rho_0 thet_0 und hlim icum]];
            %                     x0save=[x0save; x0'];
            %                     disp(x0save)
            %%
            tsave=[tsave,tt];
            if tt==1
                xsave=[xsave; x(size(x,1),:)];
                deadsave=[deadsave; x(size(x,1),12)];
            elseif tt<(maxtime-1)
                xsave=[xsave; x(size(x,1),:)];
                deadsave=[deadsave; x(size(x,1),12)];
            else
                xsave=[xsave; x(2:size(x,1),:)];
                deadsave=[deadsave; x(2:size(x,1),12)];
            end
            
            
        end
        
        savetofile(strcat(sprintf('cap_res/%d_%s',p_num,getfield(CN,regions{mm}))),tsave,xsave);
        savemulti(strcat(sprintf('rt_cap/%d_%s',p_num,getfield(CN,regions{mm}))),rtsave,parsave,x0save);
        
        %load
        
        if DECAP==1
            getthefile1=strcat(res_toPerturb,'/area_specific_',getfield(CN,regions{mm}),'.mat');
            getthefile2=strcat('rt_cap/',sprintf('%d',p_num),'_',getfield(CN,regions{mm}),'.mat');
            
            fileit=load(getthefile1);
            parfile=load(getthefile2);
            
            pop=getfield(P0,regions{mm});
            %disp(pop)
            %%
            %save parameters and initial conditions
            %             maxtime2=maxtime;
            param=parfile.parsave';%fileit.savpars;
            %                 initco=parfile.x0save;
            initco=fileit.savinit(2:size(fileit.savinit,1),:);
            
            %%
            %how to proceed:p
            %A:save initial conditions and parameters from the first window.
            ini_par=param(:,1);
            ini_init=initco(:,1);
            %D:run each; save; plot.
            
            num2=xlsread(strcat('../icu_limits/',getfield(CN,regions{mm}),'_icu.xls'));
            num3=xlsread(strcat('../hos_limits/',getfield(CN,regions{mm}),'_hos.xls'));
            
            tsave=[]; xsave=[]; rtsave=[]; parsave=[]; ddtmp=[];...
                xTestS=[]; xTestLate=[]; xTestCap=[]; xTestLateCap=[];...
                xTestS_B=[]; xTestLate_B=[]; xTestCap_B=[]; xTestLateCap_B=[];
            disp('time decap')
            disp(maxtime)
            
            %%%%TTI
            day_testEarly=7-1;
            day_testLate=21-1;%8-1+7+7
            accept=0.6;
            red_rateA=-0.02;
            red_rateB=-0.005;
            if (accept < und)
                last=accept/und;
                %     uu=undetected*[repelem(1,day_test) [1:-0.1:last] repelem(last,500)];
                reductionUND=und*[repelem(1,day_testEarly) [1:red_rateA:last] repelem(last,500)];
                UNDLate=param(15,1)*[repelem(1,day_testLate) [1:red_rateA:last] repelem(last,500)];
                reductionUND_B=und*[repelem(1,day_testEarly) [1:red_rateB:last] repelem(last,500)];
                UNDLate_B=param(15,1)*[repelem(1,day_testLate) [1:red_rateB:last] repelem(last,500)];
            else
                reductionUND=und*[repelem(1,1000)];
                UNDLate=und*[repelem(1,1000)];
                reductionUND_B=und*[repelem(1,1000)];
                UNDLate_B=und*[repelem(1,1000)];
            end
            %%%%%%%%%
            %%%%%%%%%
            for tt=1:(maxtime-1)
                if tt==1
                    tspan=0:0.2:5.2;
                elseif tt<(maxtime-1)
                    tspan=0:1;%16July
                else
                    tspan=0:6;
                end
                alpha_0=param(1,tt);
                bet_0=param(2,tt);
                delta_0=param(3,tt);
                ki_0=param(4,tt);
                r1_0=param(5,tt);
                r10_0=param(6,tt);
                r3_0=param(7,tt);
                r4_0=param(8,tt);
                r5_0=param(9,tt);
                r6_0=param(10,tt);
                r7_0=param(11,tt);
                r8_0=param(12,tt);
                rho_0=param(13,tt);
                thet_0=param(14,tt);
                und=param(15,tt);
                hlim=max(hosplim);
                icum=max(vectlim);
                
                hlimA=param(16,tt);
                icumA=param(17,tt);
                %%
                undFIX=und;
                undVAR=reductionUND(tt);
                undVARLate=UNDLate(tt);
                
                undVAR_B=reductionUND_B(tt);
                undVARLate_B=UNDLate_B(tt);
                %%
                if tt==1
                    %                         x0=ini_init;
                    x0=[getfield(P0,regions{mm}) getfield(C0,regions{mm}) repelem(0,(15-2))];
                    x0Test=[getfield(P0,regions{mm}) getfield(C0,regions{mm}) repelem(0,(16-2))];
                    x0TestLate=x0Test;
                    x0TestCap=x0Test;
                    x0TestLateCap=x0Test;
                    
                    x0Test_B=x0Test;
                    x0TestLate_B=x0Test;
                    x0TestCap_B=[getfield(P0,regions{mm}) getfield(C0,regions{mm}) repelem(0,(16-2))];
                    x0TestLateCap_B=[getfield(P0,regions{mm}) getfield(C0,regions{mm}) repelem(0,(16-2))];
                else
                    x0=xsave(tt-1,:);
                    x0Test=xTestS(tt-1,:);
                    x0TestLate=xTestLate(tt-1,:);
                    x0TestCap=xTestCap(tt-1,:);
                    x0TestLateCap=xTestLateCap(tt-1,:);
                    
                    x0Test_B=xTestS_B(tt-1,:);
                    x0TestLate_B=xTestLate_B(tt-1,:);
                    x0TestCap_B=xTestCap_B(tt-1,:);
                    x0TestLateCap_B=xTestLateCap_B(tt-1,:);
                end
                %             dbstop if naninf
                opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
                [t,x]=ode15s(@(t,x) ODEcovCAP(t,x, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    und, pop, hlim, icum),tspan,x0);%,opts);
                
                %%calc Rt
                r9=2*r3_0*r4_0/(2*r4_0+r3_0);
                m=((1-und)/(1-alpha_0));
                num=r1_0*(alpha_0/r9 + (1-alpha_0)/r3_0 + bet_0*(1-alpha_0)*(1-rho_0)*m/r4_0...
                    +bet_0*rho_0*(1-alpha_0)*m/r6_0 + ki_0*(1-alpha_0)*(1-m)/r4_0);
                frac=x(size(x,1),1)/(pop-x(size(x,1),12));
                rtsave=[rtsave; num*frac];
                parsave=[parsave;[alpha_0 bet_0 delta_0 ki_0 r1_0 r10_0 r3_0 r4_0 r5_0 r6_0 r7_0 r8_0 rho_0 thet_0 und hlim icum]];
                %%
                
                [tTest,xTest]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVAR, pop, hlim, icum),tspan,x0Test); %,opts
                [tTestL,xTestL]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVARLate, pop, hlim, icum),tspan,x0TestLate); %,opts
                
                [tTestC,xTestC]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVAR, pop, hlimA, icumA),tspan,x0TestCap); %,opts
                [tTestLC,xTestLC]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVARLate, pop, hlimA, icumA),tspan,x0TestLateCap); %,opts
                
                %%%%%%%%%%%B --> 005
                
                [tTest,xTest_B]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVAR_B, pop, hlim, icum),tspan,x0Test_B); %,opts
                [tTestL,xTestL_B]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVARLate_B, pop, hlim, icum),tspan,x0TestLate_B); %,opts
                
                [tTestC,xTestC_B]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVAR_B, pop, hlimA, icumA),tspan,x0TestCap_B); %,opts
                [tTestLC,xTestLC_B]=ode15s(@(tTest,xTest) covTEST_cap(tTest,xTest, alpha_0, bet_0, delta_0, ki_0,...
                    r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                    undFIX,undVARLate_B, pop, hlimA, icumA),tspan,x0TestLateCap_B); %,opts
                
                
                
                tsave=[tsave,tt];
                if tt==1
                    xsave=[xsave; x(size(x,1),:)];
                    xTestS=[xTestS; xTest(size(xTest,1),:)];
                    xTestLate=[xTestLate; xTestL(size(xTestL,1),:)];
                    xTestCap=[xTestCap; xTestC(size(xTestL,1),:)];
                    xTestLateCap=[xTestLateCap; xTestLC(size(xTestLC,1),:)];
                    ddtmp=[ddtmp; x(size(x,1),12)];
                    xTestS_B=[xTestS_B; xTest_B(size(xTest_B,1),:)];
                    xTestLate_B=[xTestLate_B; xTestL_B(size(xTestL_B,1),:)];
                    xTestCap_B=[xTestCap_B; xTestC_B(size(xTestL_B,1),:)];
                    xTestLateCap_B=[xTestLateCap_B; xTestLC_B(size(xTestLC_B,1),:)];
                elseif tt<(maxtime-1)
                    xsave=[xsave; x(size(x,1),:)];
                    xTestS=[xTestS; xTest(size(xTest,1),:)];
                    xTestLate=[xTestLate; xTestL(size(xTestL,1),:)];
                    xTestCap=[xTestCap; xTestC(size(xTestC,1),:)];
                    xTestLateCap=[xTestLateCap; xTestLC(size(xTestLC,1),:)];
                    ddtmp=[ddtmp; x(size(x,1),12)];
                    xTestS_B=[xTestS_B; xTest_B(size(xTest_B,1),:)];
                    xTestLate_B=[xTestLate_B; xTestL_B(size(xTestL_B,1),:)];
                    xTestCap_B=[xTestCap_B; xTestC_B(size(xTestC_B,1),:)];
                    xTestLateCap_B=[xTestLateCap_B; xTestLC_B(size(xTestLC_B,1),:)];
                else
                    xsave=[xsave; x(2:size(x,1),:)];
                    xTestS=[xTestS; xTest(2:size(xTest,1),:)];
                    xTestLate=[xTestLate; xTestL(2:size(xTestL,1),:)];
                    xTestCap=[xTestCap; xTestC(2:size(xTestC,1),:)];
                    xTestLateCap=[xTestLateCap; xTestLC(2:size(xTestLC,1),:)];
                    ddtmp=[ddtmp; x(2:size(x,1),12)];
                    xTestS_B=[xTestS_B; xTest_B(2:size(xTest_B,1),:)];
                    xTestLate_B=[xTestLate_B; xTestL_B(2:size(xTestL_B,1),:)];
                    xTestCap_B=[xTestCap_B; xTestC_B(2:size(xTestC_B,1),:)];
                    xTestLateCap_B=[xTestLateCap_B; xTestLC_B(2:size(xTestLC_B,1),:)];
                end
                %             xsave=[xsave; x(1,:)];
                %             xsave=[xsave; [x(1,:);x(size(x,1),:)]];
                
                %                     end
            end
            %         dat=[tsave',xsave];
            deaddiff=deadsave-ddtmp;
            savetofile(strcat(sprintf('dead_diff/%d_%s',p_num,getfield(CN,regions{mm}))),tsave,deaddiff);
% %             %             save(strcat(sprintf('deCap_res/%d_%s',p_num,getfield(CN,regions{mm}))),'tsave','xsave');
            savetofile(strcat(sprintf('deCap_res/%d_%s',p_num,getfield(CN,regions{mm}))),tsave,xsave);%,'.mat'
            
            savetofile(strcat(sprintf('EarlyTest_res/%d_%s_02',p_num,getfield(CN,regions{mm}))),tsave,xTestS);
            savetofile(strcat(sprintf('LateTest_res/%d_%s_02',p_num,getfield(CN,regions{mm}))),tsave,xTestLate);
            
            savetofile(strcat(sprintf('EarlyTest_res/%d_%s_cap_02',p_num,getfield(CN,regions{mm}))),tsave,xTestCap);
            savetofile(strcat(sprintf('LateTest_res/%d_%s_cap_02',p_num,getfield(CN,regions{mm}))),tsave,xTestLateCap);
            
            
            
            savetofile(strcat(sprintf('EarlyTest_res/%d_%s_005',p_num,getfield(CN,regions{mm}))),tsave,xTestS_B);
            savetofile(strcat(sprintf('LateTest_res/%d_%s_005',p_num,getfield(CN,regions{mm}))),tsave,xTestLate_B);
            
            savetofile(strcat(sprintf('EarlyTest_res/%d_%s_cap_005',p_num,getfield(CN,regions{mm}))),tsave,xTestCap_B);
            savetofile(strcat(sprintf('LateTest_res/%d_%s_cap_005',p_num,getfield(CN,regions{mm}))),tsave,xTestLateCap_B);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
    disp(getfield(CN,regions{mm}))
end
disp('SIMULATION FINISHED')

function savetofile(folder,t,x)
%     save(fullfilename,'data');
save(folder,'t','x')
end

function savemulti(folder,rtsave,parsave,x0save)
%     save(fullfilename,'data');
save(folder,'rtsave','parsave','x0save')
end


