%% this file collects the simulation results
% and generates the files to use for Fig.5
%regions={'Emilia','Piemonte','Lombardia','Liguria','Valledaosta','Marche'};
%missing=[0,0,0,1,10,2];

%% e.g.
regions={'Lombardia'};
missing=0;

for i=1:length(regions)
    dd=[];
    tp1=load(sprintf('../cap_per_results/final_params_%s.mat',regions{i}));
    for j=1:30
        dat=[];
        adjt=0;
        sim=load(sprintf('capacity_%s_fixed_%d.mat',regions{i},j));
        for k=7:6:(492-missing(i)*6)
            orgt=adjt+tp1.comb(j).toplot(1:39,k);
            orgd=tp1.comb(j).toplot(1:39,k+5);
            simd=sim.toplot(1:39,k+5);
            adjt=adjt+1;
            dat=[dat;[orgt,orgd,simd,orgd-simd]];
        end
        dd=[dd,dat];
    end
    fd=[dd(:,1),mean(dd(:,4:4:end),2),std(dd(:,4:4:end)')'];
    csvwrite(sprintf('../excess_dead_%s.csv',regions{i}),fd)
    %%this is to add in the Rscript
    disp(sprintf('%s=%d',regions{i},mean(100*(max(dd(:,4:4:120))./max(dd(:,3:4:120))))));
    
    disp(mean(max(dd(:,4:4:120))));
end

