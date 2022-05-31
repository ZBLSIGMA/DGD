function MFEA_TLS = MFEA_GHS(Tasks,pop_M,gen,rmp,p_il,reps)
%MFEA function: implementation of MFEA algorithm
    %clc    
    tic     
    pop = pop_M;
    if mod(pop,2) ~= 0
        pop = pop + 1;
    end   
    no_of_tasks=length(Tasks);
    if no_of_tasks <= 1
        error('At least 2 tasks required for MFEA');
    end
    D=zeros(1,no_of_tasks);
    for i=1:no_of_tasks
        D(i)=Tasks(i).dims;
    end
    D_multitask=max(D);
    options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','MaxIter',2);  % settings for individual learning
    M1 = ones(1,D_multitask);
    M2 = ones(1,D_multitask);
    p=0.5;
    fnceval_calls = zeros(1,reps);  
    calls_per_individual=zeros(1,pop);
    EvBestFitness = zeros(no_of_tasks*reps,gen);    % best fitness found
    TotalEvaluations=zeros(reps,gen);               % total number of task evaluations so fer
    bestobj=Inf(1,no_of_tasks);
    for rep = 1:reps
        %disp(rep)
        for i = 1 : pop
            population(i) = Chromosome();
            population(i) = initialize(population(i),D_multitask,p_il,options);
            if i < pop/2
                 population(i).skill_factor=1;
            else
                 population(i).skill_factor=2;
            end
        end
        for i = 1 : pop
            [population(i),calls_per_individual(i)] = evaluate_vec(population(i),Tasks,p_il,no_of_tasks,options);
        end
        fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
        TotalEvaluations(rep,1)=fnceval_calls(rep);
        factorial_cost=zeros(1,pop);
        for i = 1:no_of_tasks
            for j = 1:pop
                factorial_cost(j)=population(j).factorial_costs(i);
            end
            [xxx,y]=sort(factorial_cost);
            population=population(y);
            for j=1:pop
                population(j).factorial_ranks(i)=j; 
            end
            bestobj(i)=population(1).factorial_costs(i);
            EvBestFitness(i+2*(rep-1),1)=bestobj(i);
            bestInd_data(rep,i)=population(1);
        end
        
        mu = 2;     % Index of Simulated Binary Crossover (tunable)
        mum = 5;    % Index of polynomial mutation
        generation=0;
        [max_T1,max_T2,min_T1,min_T2] = cal_max_min(population);
        while generation < gen 
            population_T2=population([population.skill_factor]==2);
            generation = generation + 1;
            indorder = randperm(pop);
            if mod(generation,2) == 0
                a = 0;
            else
                a = 1;
            end
            count=1;
            for i = 1 : pop/2     
                k = 0.7 + 0.6*rand(1);
                p1 = indorder(i);
                p2 = indorder(i+(pop/2));
                child(count)=Chromosome();
                child(count+1)=Chromosome();
                child(count).survival =0;
                child(count+1).survival =0;
                u = rand(1,D_multitask);
                cf = zeros(1,D_multitask);
                cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
                cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));
                if population(p1).skill_factor == population(p2).skill_factor       % crossover      
                    if population(p1).skill_factor == 1
                        child(count) = crossover(child(count),population(p1),population(p2),cf);
                        if rand(1) > a
                            %The upper and lower limits of the unified express space are [0, 1]
                            %According to the rule of the opposite point generation in the unified representation space
                            %here is x' = 0+(1-x);
                            child(count+1).rnvec = 1 - child(count).rnvec;
                        else
                            child(count+1).rnvec = k*(max_T1+min_T1) - child(count).rnvec;
                        end
                    elseif population(p1).skill_factor == 2
                        child(count) = crossover(child(count),population(p1),population(p2),cf);
                        if rand(1) > a
                            child(count+1).rnvec = 1 - child(count).rnvec;
                        else
                            child(count+1).rnvec = k*(max_T2+min_T2) - child(count).rnvec;
                        end
                    end  
                    child(count)=mutate(child(count),child(count),D_multitask,mum);
                    child(count+1)=mutate(child(count+1),child(count+1),D_multitask,mum);
                    child(count).skill_factor = population(p1).skill_factor;
                    child(count+1).skill_factor = population(p1).skill_factor;
                elseif rand(1) < rmp
                    if rand(1) > p
                        if population(p1).skill_factor == 1 && population(p2).skill_factor == 2
                            tmp = population(p1);
                            tmp.rnvec = population(p1).rnvec .* M1;
                            tmp.rnvec(tmp.rnvec>1) = 1;
                            tmp.rnvec(tmp.rnvec<0) = 0;
                            child(count) = crossover(child(count),tmp,population(p2),cf);
                            if rand(1) > a
                                child(count+1).rnvec = 1 - child(count).rnvec;
                            else
                                child(count+1).rnvec = k*(max_T2+min_T2) - child(count).rnvec;
                            end
                        elseif population(p1).skill_factor == 2 && population(p2).skill_factor == 1
                            tmp = population(p2);
                            tmp.rnvec = population(p2).rnvec .* M1;
                            tmp.rnvec(tmp.rnvec>1) = 1;
                            tmp.rnvec(tmp.rnvec<0) = 0;
                            child(count) = crossover(child(count),population(p1),tmp,cf);
                            if rand(1) > a
                                child(count+1).rnvec = 1 - child(count).rnvec;
                            else
                                child(count+1).rnvec = k*(max_T2+min_T2) - child(count).rnvec;
                            end
                        end
                    else
                        if population(p1).skill_factor == 1 && population(p2).skill_factor == 2
                            tmp = population(p2);
                            tmp.rnvec = population(p2).rnvec .* M2;
                            tmp.rnvec(tmp.rnvec>1) = 1;
                            tmp.rnvec(tmp.rnvec<0) = 0;
                            child(count) = crossover(child(count),population(p1),tmp,cf);
                            if rand(1) > a
                                %The upper and lower limits of the unified express space are [0, 1]
                                %According to the rule of the opposite point generation in the unified representation space
                                %here is x' = 0+(1-x);
                                child(count+1).rnvec = 1 - child(count).rnvec;
                            else
                                child(count+1).rnvec = k*(max_T1+min_T1) - child(count).rnvec;
                            end          
                        elseif population(p1).skill_factor == 2 && population(p2).skill_factor == 1
                            tmp = population(p1);
                            tmp.rnvec = population(p1).rnvec .* M2;
                            tmp.rnvec(tmp.rnvec>1) = 1;
                            tmp.rnvec(tmp.rnvec<0) = 0;
                            child(count) = crossover(child(count),tmp,population(p2),cf);
                            if rand(1) > a
                                %The upper and lower limits of the unified express space are [0, 1]
                                %According to the rule of the opposite point generation in the unified representation space
                                %here is x' = 0+(1-x);
                                child(count+1).rnvec = 1 - child(count).rnvec;
                            else
                                child(count+1).rnvec = k*(max_T1+min_T1) - child(count).rnvec;
                            end
                        end
                    end
                    
                    child(count).skill_factor=round(rand(1))+1;
                    child(count+1).skill_factor=round(rand(1))+1;
                else
                    child(count)=mutate(child(count),population(p1),D_multitask,mum);
                    child(count+1)=mutate(child(count+1),population(p2),D_multitask,mum);
                    child(count).skill_factor = population(p1).skill_factor;
                    child(count+1).skill_factor = population(p2).skill_factor;
                end
                child(count).rnvec(child(count).rnvec>1)=1;
                child(count).rnvec(child(count).rnvec<0)=0;
                child(count+1).rnvec(child(count+1).rnvec>1)=1;
                child(count+1).rnvec(child(count+1).rnvec<0)=0;
                count=count+2;
            end        
            for i = 1 : pop            
                [child(i),calls_per_individual(i)] = evaluate_vec(child(i),Tasks,p_il,no_of_tasks,options);           
            end    
            
            fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
            TotalEvaluations(rep,generation)=fnceval_calls(rep);

            intpopulation(1:pop)=population;
            intpopulation(pop+1:2*pop)=child;
            factorial_cost=zeros(1,2*pop);
            for i = 1:no_of_tasks
                for j = 1:2*pop
                    factorial_cost(j)=intpopulation(j).factorial_costs(i);
                end
                [xxx,y]=sort(factorial_cost);
                intpopulation=intpopulation(y);
                for j=1:2*pop
                    intpopulation(j).factorial_ranks(i)=j;
                end
                if intpopulation(1).factorial_costs(i)<=bestobj(i)
                    bestobj(i)=intpopulation(1).factorial_costs(i);
                    bestInd_data(rep,i)=intpopulation(1);
                end
                EvBestFitness(i+2*(rep-1),generation)=bestobj(i);            
            end
            for i=1:2*pop
                [xxx,yyy]=min(intpopulation(i).factorial_ranks);
                intpopulation(i).skill_factor=yyy;
                intpopulation(i).scalar_fitness=1/xxx;
            end   
            [xxx,y]=sort(-[intpopulation.scalar_fitness]);
            intpopulation=intpopulation(y);
            population=intpopulation(1:pop);            
            [max_T1,max_T2,min_T1,min_T2] = cal_max_min(population);
            [M1,M2] = domain_ad(population,D);
            %disp(['MFEA Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj)]);         
        end 
    end
    MFEA_TLS.wall_clock_time=toc;
    MFEA_TLS.EvBestFitness=EvBestFitness;
    MFEA_TLS.bestInd_data=bestInd_data;
    MFEA_TLS.TotalEvaluations=TotalEvaluations;
end