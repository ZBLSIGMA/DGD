function [ rnvec1, rnvec2 ] = transfer( p1, p2, dim)

    p1_rnvec=	[p1(1:end).rnvec];
    p2_rnvec=	[p2(1:end).rnvec];
    p1_rnvec = vec2mat(p1_rnvec , dim);
    p2_rnvec = vec2mat(p2_rnvec , dim);
    min_point_p1 = min(p1_rnvec);
    max_point_p1 = max(p1_rnvec);
    min_point_p2 = min(p2_rnvec);
    max_point_p2 = max(p2_rnvec);
    min_all_pop = [min_point_p1;min_point_p2];
    max_all_pop = [max_point_p1;max_point_p2];
    min_all_pop = min(min_all_pop);
    max_all_pop = max(max_all_pop);
    F = 0.6;
    rand_ind = randi([1,length(p1)]);
    rand_ind_2 = randi([1,length(p1)]);
    best_ind_p1 = randi([1,length(p1)]);
    if rand(1)<1
        rnvec1 = p1(rand_ind).rnvec + F.*(p1(1).rnvec - p1(rand_ind_2).rnvec);
    else
        rnvec1 = p1(rand_ind_2).rnvec + F.*(p1(best_ind_p1).rnvec - p1(rand_ind).rnvec);
    end
    rnvec2 = max_all_pop - (rnvec1 - min_all_pop);
    rnvec1(rnvec1>1)=1;
    rnvec1(rnvec1<0)=0;
    rnvec2(rnvec2>1)=1;
    rnvec2(rnvec2<0)=0;
end

