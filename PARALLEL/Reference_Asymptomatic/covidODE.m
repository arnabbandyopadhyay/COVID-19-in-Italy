files=dir('und_results/area_specific_Italy.mat');

tspan=0:1500;
end1=109;

for i=1:length(files)
    cn=split(files(i).name,'_');
    cd('und_results/');
    k=load(files(i).name);
    cd('/home/msi18/Desktop/covid-19/COVID/undetected/')
    x0=k.savinit(2:15,end1-1);
    NI=k.savinit(2,1);
    params=num2cell(k.savpars(:,end1));
    [alpha,bet,delta,ki,r1,r10,r3,r4,r5,r6,r7,r8,rho,thet,und]=params{:};
    rt=k.savr0(end1);
%     und=0.7;
    %ki=bet;
    [t,x]=ode45(@(t,x) covid19ODE(t,x,alpha,bet,delta,ki,r1,r10,r3,r4,r5,r6,r7,r8,rho,thet,und,rt,NI),tspan,x0);
    dat=[t,x];
    save(strcat('extension_',cn{3}),'t','x'); %,'.mat'
    
end
plot(t(:,1),x(:,1))
figure
plot(t(2:length(t),1),diff(x(:,1)))
function dydt = covid19ODE(t,x,alpha,bet,delta,ki,r1,r3,r4,r5,r6,r7,r8,r10,rho,thet,und,NI)

N = 14;
dydt = zeros(N,1);


dydt(1) =  -r1*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)))*x(1)/NI; % sus
dydt(2) =  r1*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)))*x(1)/NI -(r3/(r3*5.2 -1))*x(2); %exposed
dydt(3) =  (1-alpha)*(r3/(r3*5.2 -1))*x(2) - r3*x(3); %car
dydt(4) = alpha*(r3/(r3*5.2 -1))*x(2) - x(4)*(2*r3*r4/(2*r4+r3)); % car recovered
dydt(5) = ((1-und)/(1-alpha))*rho*r3*x(3)-r6*x(5); % inf
dydt(6) = ((1-und)/(1-alpha))*(1-rho)*r3*x(3)-r4*x(6); % inf recovered
dydt(7) = thet*r6*x(5) - r7*x(7); % hos
dydt(8) = (1-thet)*r6*x(5) - x(8)*r5; % hos recovered
dydt(9) = delta*r7*x(7) - r10*x(9); % icu
dydt(10) = (1-delta)*r7*x(7) - x(10)*r8; % icu recovered
dydt(11) = x(8)*r5 + x(10)*r8; % rec
dydt(12) = r10*x(9); % dead
dydt(13) = x(6)*r4 + x(4)*(2*r3*r4/(2*r4+r3))+r4*x(14); % rec residual
dydt(14) = (1-((1-und)/(1-alpha)))*r3*x(3)-r4*x(14); % inf Undetected
end
