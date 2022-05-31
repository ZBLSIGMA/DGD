function [M1,M2] = domain_ad(population,D)
    population_T1=population([population.skill_factor]==1);
    population_T2=population([population.skill_factor]==2);
    T1 = [];
    T2 = [];
    
    N1 = unidrnd(2);
    for i=1:N1
        T1 = [T1;population_T1(i).rnvec];
    end
    for i=1:N1
        T2 = [T2;population_T2(i).rnvec];
    end
    mean_T1 = mean(T1);
    mean_T2 = mean(T2);
%     M1 = mean_T2 ./ mean_T1;
%     M2 = mean_T1 ./ mean_T2;
    M1 = (mean_T2+0.0000000000000000001) ./ (mean_T1+0.0000000000000000001);
    M2 = (mean_T1+0.0000000000000000001) ./ (mean_T2+0.0000000000000000001);
end 