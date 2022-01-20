function ar=exp_loop(datatoread,deffile,loopfile,cityname,missing,undetected,home_dir)

shift=5.2;
datapoints=7;
parstosave=17;

mkdir(strcat(cityname))
mkdir(strcat(strcat(cityname),'/Models'))
mkdir(strcat(strcat(cityname),'/Data'))
copyfile(strcat('Models/',strcat(deffile,'.def')), strcat(strcat(cityname),'/Models'))
copyfile(strcat('Models/',loopfile), strcat(strcat(cityname),'/Models'))
copyfile('xlwrite.m',strcat(cityname))
copyfile('Copy_of_func_replace_string.m',strcat(cityname))
copyfile('R0calc.m',strcat(cityname))
%take the raw dataset
[num,txt,raw]=xlsread(datatoread);

copyfile(strcat('icu_limits/',strcat(cityname,'_icu.xls')),strcat(cityname));
copyfile(strcat('hos_limits/',strcat(cityname,'_hos.xls')),strcat(cityname));

cd(strcat(cityname))

num2=xlsread(strcat(cityname,'_icu.xls'));
num3=xlsread(strcat(cityname,'_hos.xls'));

if missing>0
    %         vectlim=num2((missing+1):length(num2),2);
    %         hosplim=num3((missing+1):length(num3),2);
    vectlim=num2((missing+4):length(num2),2);
    hosplim=num3((missing+4):length(num3),2);
else
    %         vectlim=num2(1:length(num2),2);
    %         hosplim=num3(1:length(num3),2);
    vectlim=num2(1+3:length(num2),2);
    hosplim=num3(1+3:length(num3),2);
end

seq=[3:length(raw)-datapoints+1];
%     seq=[3:55]; % for testing
disp(cityname)
disp(size(seq))
end
