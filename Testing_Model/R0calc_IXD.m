function R0 = R0calc_IXD(params,label,uu,frac)

r1=params(find(strcmp(label,'r1')));alpha=params(find(strcmp(label,'alpha')));
r3=params(find(strcmp(label,'r3')));rho=params(find(strcmp(label,'rho')));
r4=params(find(strcmp(label,'r4')));bet=params(find(strcmp(label,'bet')));
r6=params(find(strcmp(label,'r6')));und=params(find(strcmp(label,'und')));

m=((1-uu)/(1-alpha));

gamma_p = (und-alpha)/(1-alpha);

delta_p = (uu-und)/(1-alpha);

ki=params(find(strcmp(label,'ki')));

r9=2*r3*r4/(2*r4+r3);
    
num=r1*(alpha/r9 + (1-alpha)/r3 + bet*(1-alpha)*(1-rho)*m/r4 +bet*rho*(1-alpha)*m/r6 + ...
    bet*(1-alpha)*delta_p/r4 + ki*(1-alpha)*gamma_p/r4);

R0=num*frac;
    
end
