function dydt = ODEcovCAP(t,x, alpha, bet, delta, ki,...
            r1, r10, r3, r4, r5, r6, r7, r8, rho, thet,...
            und, pop, hlim,icum)

        
N=15;
dydt = zeros(N,1);
%%
dydt(1) = -r1*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)+x(15)))*x(1)/pop;  % sus _a1 x1
dydt(2) = r1*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)+x(15)))*x(1)/pop -(r3/(r3*5.2 -1))*x(2);  % exp x2
dydt(3) = (1-alpha)*(r3/(r3*5.2 -1))*x(2) - r3*x(3);  % CI_a1 x3
dydt(4) = alpha*(r3/(r3*5.2 -1))*x(2) - x(4)*(2*r3*r4/(2*r4+r3));  % CR x4
dydt(5) = ((1-und)/(1-alpha))*rho*r3*x(3)-r6*x(5)*(1-(exp(x(7)+x(8)-hlim)^10/(1+exp(x(7)+x(8)-hlim)^10))^10) - r6*x(5)*((exp(x(7)+x(8)-hlim)^10/(1+exp(x(7)+x(8)-hlim)^10))^10);  % IH x5
dydt(6) = ((1-und)/(1-alpha))*(1-rho)*r3*x(3)-r4*x(6);  % IR x6
dydt(7) = thet*r6*x(5)*(1-(exp(x(7)+x(8)-hlim)^10/(1+exp(x(7)+x(8)-hlim)^10))^10) - r7*x(7);  % HU x7
dydt(8) = (1-thet)*r6*x(5)*(1-(exp(x(7)+x(8)-hlim)^10/(1+exp(x(7)+x(8)-hlim)^10))^10) - x(8)*r5;  % HR x8
dydt(9) = delta*r7*x(7)*(1-(exp(x(9)+x(10)-icum)^10/(1+exp(x(9)+x(10)-icum)^10))^10) - r10*x(9);  %  UD x9
dydt(10) = (1-delta)*r7*x(7)*(1-(exp(x(9)+x(10)-icum)^10/(1+exp(x(9)+x(10)-icum)^10))^10) - x(10)*r8;  % UR x10
dydt(11) = x(6)*r4 + x(8)*r5 + x(10)*r8;  %   RZ x11
dydt(12) = r10*x(9) + r7*x(7)*(exp(x(9)+x(10)-icum)^10/(1+exp(x(9)+x(10)-icum)^10))^10 + r7*x(15);  %  dead x12
dydt(13) = x(4)*(2*r3*r4/(2*r4+r3))+r4*x(14);  %  RX x13
dydt(14) = (1-((1-und)/(1-alpha)))*r3*x(3)-r4*x(14);  %  IX x14
dydt(15) = r6*x(5)*((exp(x(7)+x(8)-hlim)^10/(1+exp(x(7)+x(8)-hlim)^10))^10) - r7*x(15); %infDead x15



end
