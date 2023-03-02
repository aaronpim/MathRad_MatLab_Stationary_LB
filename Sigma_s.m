function [Sctr] = Sigma_s(x,s1,s2)
    [X,Y] = meshgrid(x,x);
    R = sqrt(X.^2+Y.^2);
    
    s_dot_s = -5:5;
    N = length(s_dot_s);
    
    prob_vec_inclusions = ones(N,1); 
    prob_vec_inclusions = prob_vec_inclusions./sum(prob_vec_inclusions);
    prob_vec_exclusions = [1;zeros(N-1,1)]; 
    prob_vec_exclusions = prob_vec_exclusions./sum(prob_vec_exclusions);
    
    Sctr = (R>0.5).*dot(prob_vec_exclusions, s_dot_s==dot(s1,s2))...
         +(R<=0.5).*dot(prob_vec_inclusions, s_dot_s==dot(s1,s2));
end

