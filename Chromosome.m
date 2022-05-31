classdef Chromosome    
    properties
        rnvec; % (genotype)--> decode to find design variables --> (phenotype) 
        factorial_costs;
        factorial_ranks;
        scalar_fitness;
        skill_factor;
        survival=0;
        p_il;
        options;
    end    
    methods        
        function object = initialize(object,D,p_il,options)            
            object.rnvec = rand(1,D);
            object.p_il = p_il;
            object.options = options;
        end
        
        function [object,calls] = evaluate_vec(object,Tasks,p_il,no_of_tasks,options)
            calls = 0;
            if object.skill_factor == 0
                calls=0;
                for i = 1:no_of_tasks
                    [object.factorial_costs(i),xxx,funcCount]=fnceval(Tasks(i),object.rnvec,p_il,options);
                    calls = calls + funcCount;
                end
            else
                object.factorial_costs(1:no_of_tasks)=inf;
                %object.factorial_costs(no_of_tasks)=inf;
                for i = 1:no_of_tasks
                    if object.skill_factor == i
                        [object.factorial_costs(object.skill_factor),object.rnvec,funcCount]=fnceval(Tasks(object.skill_factor),object.rnvec,p_il,options);
                        calls = funcCount;
                        break;
                    end
                end
            end
        end
        
        function [object,calls] = evaluate_SOO(object,Task,p_il,options)   
            [object.factorial_costs,object.rnvec,funcCount]=fnceval(Task,object.rnvec,p_il,options);
            calls = funcCount;
        end
        
        % SBX
        function object=crossover(object,p1,p2,cf)
            object.rnvec=0.5*((1+cf).*p1.rnvec + (1-cf).*p2.rnvec);
            object.rnvec(object.rnvec>1)=1;
            object.rnvec(object.rnvec<0)=0;
        end

        function object=crossover11(object,p1,p2,cf,QWE,L,sigma)             
            p1.rnvec=p1.rnvec-QWE(1,:).*sigma/L;
             p2.rnvec=p2.rnvec-QWE(2,:).*sigma/L;
            object.rnvec=0.5*((1+cf).*p1.rnvec + (1-cf).*p2.rnvec);
            object.rnvec(object.rnvec>1)=1;
            object.rnvec(object.rnvec<0)=0;
        end

        function object=mutate1(object,p,dim,qq,L,sigma)          
            object.rnvec = p.rnvec-qq.*sigma/L;          
        end    
        
        % polynomial mutation
        function object=mutate(object,p,dim,mum)
            rnvec_temp=p.rnvec;
            for i=1:dim
                if rand(1)<1/dim
                    u=rand(1);
                    if u <= 0.5
                        del=(2*u)^(1/(1+mum)) - 1;
                        rnvec_temp(i)=p.rnvec(i) + del*(p.rnvec(i));
                    else
                        del= 1 - (2*(1-u))^(1/(1+mum));
                        rnvec_temp(i)=p.rnvec(i) + del*(1-p.rnvec(i));
                    end
                end
            end  
            object.rnvec = rnvec_temp;          
        end    
    end
end