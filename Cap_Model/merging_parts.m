
regions={'Lombardia'};

for i=1:length(regions)
    disp(i)
    
        
    tp1=load(sprintf('area_specific_%s_part1.mat',regions{i}));
    tp2=load(sprintf('area_specific_%s_part2.mat',regions{i}));
    tp3=load(sprintf('area_specific_%s_part3.mat',regions{i}));
    tp4=load(sprintf('area_specific_%s_part4.mat',regions{i}));
    tp5=load(sprintf('area_specific_%s_part5.mat',regions{i}));
    tp6=load(sprintf('area_specific_%s_part6.mat',regions{i}));
    
    comb=[tp1.TP,tp2.TP,tp3.TP,tp4.TP,tp5.TP,tp6.TP];
    save([sprintf('final_params_%s.mat',regions{i})],'comb');
    
        
    
end
