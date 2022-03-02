function dydt = covid19ODE(t,x,alpha,bet,delta,ki,r1,r10,r3,r4,r5,r6,r7,r8,rho,thet,und,rt,NI)

N = 14;
dydt = zeros(N,1);
m=((1-und)/(1-alpha));
r9=2*r3*r4/(2*r4+r3);

A=(alpha/r9 + (1-alpha)/r3 + bet*(1-alpha)*(1-rho)*m/r4 +bet*rho*(1-alpha)*m/r6 + ki*(1-alpha)*(1-m)/r4);
ffrac = rt/(r1*A);

r1=1/(A);

dydt(1) =  -(r1)*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6))); % sus
dydt(2) =  (r1)*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6))) -(r3/(r3*5.2 -1))*x(2); %exposed
% dydt(1) = -rt*((NI-x(12))/(A*x(1)))*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)))*x(1)/(NI-x(12)); % sus
% dydt(2) = rt*((NI-x(12))/(A*x(1)))*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)))*x(1)/(NI-x(12)) -(r3/(r3*5.2 -1))*x(2); %exposed
dydt(3) = (1-alpha)*(r3/(r3*5.2 -1))*x(2) - r3*x(3); %car
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
