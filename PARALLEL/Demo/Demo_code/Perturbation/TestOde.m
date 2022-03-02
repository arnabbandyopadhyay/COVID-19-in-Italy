function dydt = covid19ODE(t,x, alpha, bet, delta, ki,...
    r1, r10, r3, r4, r5, r6, r7, r8, rho, thet,...
    undFIX, undVAR, pop)


%    1 {'alpha'}    2{'bet'}    3-6{'d'}x4    7-10{'delta'}x4
%     11{'ki'}    12-15{'r1'}x4    16{'r3'}    17{'r4'}    18{'r5'}
%    19{'r6'}    20{'r7'}    21{'r8'}    22-25{'rho'}x4    26-29{'thet'}x4
%     30{'und'}   31{'zcm11'}    32{'zcm12'}    33{'zcm13'}
%     34{'zcm14'}    35{'zcm21'}    36{'zcm22'}    37{'zcm23'}
%     38{'zcm24'}    39{'zcm31'}    40{'zcm32'}    41{'zcm33'}
%     42{'zcm34'}    43{'zcm41'}    44{'zcm42'}
%     45{'zcm43'}    46{'zcm44'}

N=15;
dydt = zeros(N,1);
%%
dydt(1) = -r1*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)+x(15)))*x(1)/pop;  % sus _a1 x1
dydt(2) = r1*(x(3)+x(4)+ki*x(14)+bet*(x(5)+x(6)+x(15)))*x(1)/pop -(r3/(r3*5.2 -1))*x(2);  % exp x2
dydt(3) = (1-alpha)*(r3/(r3*5.2 -1))*x(2) - r3*x(3);  % CI_a1 x3
dydt(4) = alpha*(r3/(r3*5.2 -1))*x(2) - x(4)*(2*r3*r4/(2*r4+r3));  % CR x4
dydt(5) = ((1-undFIX)/(1-alpha))*rho*r3*x(3)-r6*x(5);  % IH x5
dydt(6) = ((1-undFIX)/(1-alpha))*(1-rho)*r3*x(3)-r4*x(6);  % IR x6
dydt(7) = thet*r6*x(5) - r7*x(7);  % HU x7
dydt(8) = (1-thet)*r6*x(5) - x(8)*r5;  % HR x8
dydt(9) = delta*r7*x(7) - r10*x(9);  %  UD x9
dydt(10) = (1-delta)*r7*x(7) - x(10)*r8;  % UR x10
dydt(11) = r4*x(15)+x(6)*r4 + x(8)*r5 + x(10)*r8;  %   RZ x11
dydt(12) = r10*x(9);  %  dead x12
dydt(13) = x(4)*(2*r3*r4/(2*r4+r3))+r4*x(14) + r4*x(14);  %  RX x13
dydt(14) = ((undVAR-alpha)/(1-alpha))*r3*x(3)-r4*x(14);  %  IX x14
dydt(15) = ((undFIX - undVAR)/(1-alpha))*r3*x(3) -r4*x(15); %Ixd x15


end

%
%
%
%
% "-r1*(x3+x4+ki*x14+bet*(x5+x6+x15))*x1/60359546"     // sus x1
% "r1*(x3+x4+ki*x14+bet*(x5+x6+x15))*x1/60359546 -(r3/(r3*5.2 -1))*x2"     // exposed x2
% "(1-alpha)*(r3/(r3*5.2 -1))*x2 - r3*x3" //car x3
% "alpha*(r3/(r3*5.2 -1))*x2 - x4*(2*r3*r4/(2*r4+r3))" //car2 x4
% "((1-9.357000e-01)/(1-alpha))*rho*r3*x3-r6*x5" //inf x5 use r6
% "((1-9.357000e-01)/(1-alpha))*(1-rho)*r3*x3-r4*x6" //inf2 x6
% "thet*r6*x5 - r7*x7" //hos x7 use r7
% "(1-thet)*r6*x5 - x8*r5" //hos2 x8
% "delta*r7*x7 - r10*x9" // Icu x9
% "(1-delta)*r7*x7 - x10*r8" //icu2 x10 use r8
% "x8*r5 + x10*r8"  // rec x11
% "r10*x9"  //dead x12 use d
% "x6*r4 + x4*(2*r3*r4/(2*r4+r3))+r4*x14"  //rec residual
% "((und-alpha)/(1-alpha))*r3*x3-r4*x14" //infUND x14 gamma*Ci
% ((undVAR-alpha)/(1-alpha))*r3*x(3)-r4*x(14)
% "((9.357000e-01 - und)/(1-alpha))*r3*x3 -r4*x15"     // Ixd x15
