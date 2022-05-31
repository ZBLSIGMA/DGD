function data_MFEA = MFEAAKT(Tasks,pop,gen,rmp,p_il,reps)
%MFEA function: implementation of MFEA algorithm
%clc
tic
if mod(pop,2) ~= 0
    pop = pop + 1;
end
no_of_tasks=length(Tasks);
%no_of_tasks = 2;
if no_of_tasks <= 1
    error('At least 2 tasks required for MFEA');
end
D=zeros(1,no_of_tasks);
for i=1:no_of_tasks
    D(i)=Tasks(i).dims;
end
D_multitask=max(D);
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','MaxIter',2);  % settings for individual learning
fnceval_calls = zeros(1,reps);
calls_per_individual=zeros(1,pop);
EvBestFitness = zeros(no_of_tasks*reps,gen);    % best fitness found
TotalEvaluations=zeros(reps,gen);               % total number of task evaluations so fer
bestobj=Inf(1,no_of_tasks);
for rep = 1:reps
    skill_factor = 1;
    subpop=pop/no_of_tasks;
    count =1;
    for t = 1:no_of_tasks
        for i = 1 : pop/no_of_tasks
            population(count) = Chromosome();
            population(count) = initialize(population(count),D_multitask);
            population(count).skill_factor=skill_factor;
            count = count+1;
        end
        skill_factor = skill_factor+1;
    end
    for i = 1 : pop
        child(i)=Chromosome();
        %[population(i)] = evaluate(population(i),Tasks);
        [population(i)] = evaluate_vec(population(i),Tasks,p_il,no_of_tasks,options);
    end
    
    fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
    TotalEvaluations(rep,1)=fnceval_calls(rep);
    
    for i = 1:no_of_tasks
        temp_pop = population((i-1)*subpop+1 : i*subpop);
        factorial_cost = [temp_pop.factorial_costs];
        %disp(factorial_cost);
        factorial_cost_single = factorial_cost(i:2:end);
        %disp(factorial_cost);
        %input('s')
        [xxx,y]=sort(factorial_cost_single);
        %disp(y);
        temp_pop=temp_pop(y);
        for j=1:subpop
            temp_pop(j).factorial_ranks=j;
        end
        obj_two = temp_pop(1).factorial_costs;
        obj_one = obj_two(i);
        %disp(obj_one);
        bestobj(i)=obj_one; %temp_pop(1).factorial_costs;
        EvBestFitness(i+no_of_tasks*(rep-1),1)=bestobj(i);
        bestInd_data(rep,i)=temp_pop(1);
        population((i-1)*subpop+1 : i*subpop) = temp_pop;
    end
    mu = 2;     % Index of Simulated Binary Crossover (tunable)
    mum = 20;    % Index of polynomial mutation
    generation=1;
    while generation < gen
        generation = generation + 1;
        indorder = randperm(pop);
        count=1;
        for i = 1 : pop/2
            p1 = indorder(i);
            p2 = indorder(i+(pop/2));
            u = rand(1,D_multitask);
            cf = zeros(1,D_multitask);
            cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
            cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));
            if (population(p1).skill_factor == population(p2).skill_factor) || (rand(1)<rmp)      % crossover
                pop1=population((population(p1).skill_factor-1)*subpop + 1: population(p1).skill_factor*subpop);
                pop2=population((population(p2).skill_factor-1)*subpop + 1: population(p2).skill_factor*subpop);
                if rand(1)<0.5
                    [child(count).rnvec, child(count+1).rnvec]=transfer( pop1, pop2, D_multitask);
                else
                    child(count) = crossover(child(count),population(p1),population(p2),cf);
                    child(count+1) = crossover(child(count+1),population(p2),population(p1),cf);
                    child(count)=mutate(child(count),child(count),D_multitask,mum);
                    child(count+1)=mutate(child(count+1),child(count+1),D_multitask,mum);
                end
                sf1=1+round(rand(1));
                sf2=1+round(rand(1));
                if sf1 == 1 % skill factor selection
                    child(count).skill_factor=population(p1).skill_factor;
                else
                    child(count).skill_factor=population(p2).skill_factor;
                end
                
                if sf2 == 1
                    child(count+1).skill_factor=population(p1).skill_factor;
                else
                    child(count+1).skill_factor=population(p2).skill_factor;
                end
                
            else
                child(count)=mutate(child(count),population(p1),2,mum);
                child(count).skill_factor=population(p1).skill_factor;
                child(count+1)=mutate(child(count+1),population(p2),2,mum);
                child(count+1).skill_factor=population(p2).skill_factor;
            end
            count=count+2;
        end
        for i = 1 : pop
            [child(i)] = evaluate_vec(child(i),Tasks,p_il,no_of_tasks,options);
        end
        fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
        TotalEvaluations(rep,generation)=fnceval_calls(rep);
        
        intpopulation(1:pop)=population;
        intpopulation(pop+1:2*pop)=child;
        for i = 1:no_of_tasks
            temppopulation = intpopulation([intpopulation.skill_factor]==i);
            factorial_cost = [temppopulation.factorial_costs];
            [xxx,y]=sort(factorial_cost);
            temppopulation=temppopulation(y);
            for j=1:length(temppopulation)
                temppopulation(j).factorial_ranks=j;
                temppopulation(j).scalar_fitness = 1/j;
            end
            if temppopulation(1).factorial_costs<=bestobj(i)
                bestobj(i)=temppopulation(1).factorial_costs;
                bestInd_data(rep,i)=temppopulation(1);
            end
            population((i-1)*subpop + 1: i*subpop) = temppopulation(1:subpop);
            EvBestFitness(i+no_of_tasks*(rep-1),generation)=bestobj(i);
        end
        
        disp(['MFEA Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj)]);
    end
end
data_MFEA.wall_clock_time=toc;
data_MFEA.EvBestFitness=EvBestFitness;
data_MFEA.bestInd_data=bestInd_data;
data_MFEA.TotalEvaluations=TotalEvaluations;
end