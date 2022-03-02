% close all
clear all
clc

num_pert=100;
var_toApply=0.1;

mkdir('Asym_pert_res')
mkdir('Asym_rt_res')


mainfol=pwd;

regions = [{'italy'}];


P0.abruzzo=1311580;P0.campania=5801692;P0.emilia=4459477;P0.friuli=1215220;
P0.lazio=5879082;P0.liguria=1550640;P0.lombardia=10060574;P0.italy=60359546;
P0.marche=1525271;P0.molise=305617;P0.piemonte=4356406;P0.sardegna=1639591;
P0.sicilia=4999891;P0.toscana=3729641;P0.valledaosta=125666;P0.veneto=4905854;
P0.calabria=1947131;P0.puglia=4029053;P0.umbria=882015;P0.basilicata=562869;
P0.bolzano=520891;P0.trento=538223;

%%%initial conditions: only asymptomatic model
uu=0.4;
C0.abruzzo=1*1/(1-uu);C0.campania=3*1/(1-uu);C0.emilia=18*1/(1-uu);C0.friuli=6*1/(1-uu);
C0.italy=221*1/(1-uu);C0.lazio=3*1/(1-uu);C0.liguria=1*1/(1-uu);C0.lombardia=225*1/(1-uu);
C0.marche=1*1/(1-uu);C0.molise=3*1/(1-uu);C0.piemonte=2*1/(1-uu);C0.sardegna=1*1/(1-uu);
C0.sicilia=1*1/(1-uu);C0.toscana=2*1/(1-uu);C0.valledaosta=2*1/(1-uu);C0.veneto=32*1/(1-uu);
C0.calabria=1*1/(1-uu);C0.puglia=1*1/(1-uu);C0.umbria=2*1/(1-uu);C0.basilicata=1*1/(1-uu);
C0.bolzano=1*1/(1-uu);C0.trento=4*1/(1-uu);

CN.abruzzo='Abruzzo'; CN.campania='Campania'; CN.emilia='Emilia'; CN.friuli='Friuli';
CN.italy='Italy'; CN.lazio='Lazio'; CN.liguria='Liguria'; CN.lombardia='Lombardia';
CN.marche='Marche'; CN.molise='Molise'; CN.piemonte='Piemonte'; CN.sardegna='Sardegna';
CN.sicilia='Sicilia'; CN.toscana='Toscana'; CN.valledaosta='Valledaosta'; CN.veneto='Veneto';
CN.calabria='Calabria'; CN.puglia='Puglia'; CN.umbria='Umbria'; CN.basilicata='Basilicata';
CN.bolzano='Bolzano'; CN.trento='Trento';


parfor mm=1:size(regions,2)%numel(fieldnames(CN))
    
    %%
    %load
    getthefile=strcat('Asymp_results/area_specific_',getfield(CN,regions{mm}),'.mat');
    fileit=load(getthefile);
    
    disp(getfield(CN,regions{mm}))
    pop=getfield(P0,regions{mm});
    %disp(pop)
    %%
    %save parameters and initial conditions
    CHOOSEINI=1;
    param=fileit.savpars;
    initco=fileit.savinit(2:size(fileit.savinit,1),:);
    maxtime=(size(param,2)-1);
    disp(maxtime)
    %%
    %how to proceed:
    %A:save initial conditions and parameters from the first window.
    ini_par=param(:,CHOOSEINI);
    ini_init=initco(:,1);
    
    r1r=[(0.001),5];%5
    r3r=[(1/4.2),(2/5.2)];%7
    r4r=[(1/14),(1/7)];%8
    
    r5r=[(1/16),(1/5)];%9
    r6r=[(1/7),(0.9)];%10
    r7r=[(1/3.5),(1)];%11
    r8r=[(1/16),(1/3)];%12
    r10r=[(0.1),(0.9)];%6
    
    alphar=[(0.01),(0.4)];%1
    
    betr=[(0.05),(0.1)];%2
    deltar=[(0.1),(0.9)];%3
    
    rhor=[(0.01),(0.9)];%13
    thetr=[(0.01),(0.7)];%14
    
    kir=[(10^-8),(1)];%4
    undr=[(0.002),(0.99)];%15
    
    alpha_0=ini_par(1);
    bet_0=ini_par(2);
    delta_0=ini_par(3);
    ki_0=ini_par(4);
    r1_0=ini_par(5);
    r10_0=ini_par(6);
    r3_0=ini_par(7);
    r4_0=ini_par(8);
    r5_0=ini_par(9);
    r6_0=ini_par(10);
    r7_0=ini_par(11);
    r8_0=ini_par(12);
    rho_0=ini_par(13);
    thet_0=ini_par(14);
    und=ini_par(15);
    
    rhos=param(13,CHOOSEINI:maxtime);
    thets=param(14,CHOOSEINI:maxtime);
    r1s=param(5,CHOOSEINI:maxtime);
    r10s=param(6,CHOOSEINI:maxtime);
    deltas=param(3,CHOOSEINI:maxtime);
    
    %%parameter bounds
    lbper=[alphar(1),betr(1),deltar(1),kir(1),r1r(1),r10r(1),r3r(1),...
        r4r(1),r5r(1),r6r(1),r7r(1),r8r(1),rhor(1),thetr(1),undr(1)];
    ubper=[alphar(2),betr(2),deltar(2),kir(2),r1r(2),r10r(2),r3r(2),...
        r4r(2),r5r(2),r6r(2),r7r(2),r8r(2),rhor(2),thetr(2),undr(2)];
    
    lbper2=[rhor(1),thetr(1),r10r(1),deltar(1),r1r(1)];
    ubper2=[rhor(2),thetr(2),r10r(2),deltar(2),r1r(2)];
    
    parchar=["rho","thet","r10","delta","r1"];
    
    
    %%
    %B:apply +/-10% to the fractions (und excluded) WITHIN their range.
    rhoperF=[]; thetperF=[]; r10perF=[]; deltaperF=[]; r1perF=[];
    for outind=CHOOSEINI:maxtime
        disp(outind)
        
        k0=rhos(outind)*thets(outind)*r10s(outind)*deltas(outind)*r1s(outind);
        base_par=[rhos(outind),thets(outind),r10s(outind),deltas(outind),r1s(outind)];
        %%%
        ini_par=param;%(:,2); ???MS
        
        plbper=base_par-base_par.*var_toApply;
        pubper=base_par+base_par.*var_toApply;
        
        for parlist=1:size(plbper,1)
            if plbper(parlist)<lbper2(parlist)
                plbper(parlist)=lbper2(parlist);
            end
            if pubper(parlist)>ubper2(parlist)
                pubper(parlist)=ubper2(parlist);
            end
        end
        
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
    
    %D:run each; save.
    
    %%%%%%%%%
    
    for p_num=1:(sizelt+1)
        tsave=[]; xsave=[]; rtsave=[]; parsave=[]; xTestS=[]; xTestLate=[];
        rtearly=[]; rtlate=[];
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
            
            if tt==1
                x0=[getfield(P0,regions{mm}) getfield(C0,regions{mm}) repelem(0,(14-2))];
            else
                x0=xsave(tt-1,:);
            end
            
            [t,x]=ode15s(@(t,x) covid19ODE(t,x, alpha_0, bet_0, delta_0, ki_0,...
                r1_0, r10_0, r3_0, r4_0, r5_0, r6_0, r7_0, r8_0, rho_0, thet_0,...
                und, pop),tspan,x0);
            %%calc Rt
            r9=2*r3_0*r4_0/(2*r4_0+r3_0);
            m=((1-und)/(1-alpha_0));
            num=r1_0*(alpha_0/r9 + (1-alpha_0)/r3_0 + bet_0*(1-alpha_0)*(1-rho_0)*m/r4_0...
                +bet_0*rho_0*(1-alpha_0)*m/r6_0 + ki_0*(1-alpha_0)*(1-m)/r4_0);
            frac=x(size(x,1),1)/(pop-x(size(x,1),12));
            rtsave=[rtsave; num*frac];
            parsave=[parsave;[alpha_0 bet_0 delta_0 ki_0 r1_0 r10_0 r3_0 r4_0 r5_0 r6_0 r7_0 r8_0 rho_0 thet_0 und]];
            %%
            %%%
            %%
            
            %%%
            tsave=[tsave,tt];
            if tt==1
                xsave=[xsave; x(size(x,1),:)];
            elseif tt<(maxtime-1)
                xsave=[xsave; x(size(x,1),:)];
            else
                xsave=[xsave; x(2:size(x,1),:)];
            end
            
        end
        %         dat=[tsave',xsave];
%         save(strcat(sprintf('Asym_pert_res/%d_%s',p_num,getfield(CN,regions{mm}))),'tsave','xsave'); %,'.mat'
        savetofile(strcat(sprintf('Asym_pert_res/%d_%s',p_num,getfield(CN,regions{mm}))),tsave,xsave);
        
%         save(strcat(sprintf('Asym_rt_res/%d_%s',p_num,getfield(CN,regions{mm}))),'rtsave','parsave'); %,'.mat'
        savetofile(strcat(sprintf('Asym_rt_res/%d_%s',p_num,getfield(CN,regions{mm}))),rtsave,parsave);
    end
    
end

disp('DONE')


function savetofile(folder,t,x)
%     save(fullfilename,'data');
save(folder,'t','x')
end

