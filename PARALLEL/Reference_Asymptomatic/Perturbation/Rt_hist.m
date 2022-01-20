rt=[];
for i=1:1001
    load((sprintf('%d_Italy.mat',i)))
    rt=[rt; t(100)];
end

hist(rt)