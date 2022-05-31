function [max_T1,max_T2,min_T1,min_T2] = cal_max_min(population)
    rnvec_T1 = [];
    rnvec_T2 = [];
    population_T1=population([population.skill_factor]==1);
    population_T2=population([population.skill_factor]==2);
    for i =1:length(population_T1)
        rnvec_T1 = [rnvec_T1;population_T1(i).rnvec];
    end
    for i =1:length(population_T2)
        rnvec_T2 = [rnvec_T2;population_T2(i).rnvec];
    end
    max_T1 = max(rnvec_T1);
    max_T2 = max(rnvec_T2);
    min_T1 = min(rnvec_T1);
    min_T2 = min(rnvec_T2);
end