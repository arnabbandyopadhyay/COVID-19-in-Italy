function R0 = R0calc(params,label,frac)

r1=params(find(strcmp(label,'r1')));alpha=params(find(strcmp(label,'alpha')));
r3=params(find(strcmp(label,'r3')));rho=params(find(strcmp(label,'rho')));
r4=params(find(strcmp(label,'r4')));bet=params(find(strcmp(label,'bet')));
r6=params(find(strcmp(label,'r6')));und=params(find(strcmp(label,'und')));
m=((1-und)/(1-alpha));%params(find(strcmp(label,'m')));
ki=params(find(strcmp(label,'ki')));

% r6=r4+params(find(strcmp(label,'r66')));
r9=2*r3*r4/(2*r4+r3);

% num=params(6) * ((1-params(15))*params(10) + params(2)*(1-params(1))*params(8) + ...
%     params(15)*params(12));
    
num=r1*(alpha/r9 + (1-alpha)/r3 + bet*(1-alpha)*(1-rho)*m/r4 +bet*rho*(1-alpha)*m/r6 + ki*(1-alpha)*(1-m)/r4);

% den=((1-params(1))*params(8) + ...
%     params(1)*params(9) ...
%     ) * ((1-params(15))*params(10) + ...
%     params(15)*params(12));
    %coeff=S0/N0 //???????
R0=num*frac;
    
end