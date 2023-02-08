function [Sctr] = Sigma_s(x,s1,s2)
[X,Y,Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);

dot_vec = [-1,-2/sqrt(6),-1/sqrt(2),-1/sqrt(3),-1/2,-1/3,0,1/3,1/2,1/sqrt(3),1/sqrt(2),2/sqrt(6),1];

%prob_vec_inclusions = (1:length(dot_vec))./sum(1:length(dot_vec));
prob_vec_inclusions = ones(length(dot_vec),1)/length(dot_vec);
prob_vec_exclusions = [zeros(1,11),0.1,0.9];

Sctr = (R>0.5).*dot(prob_vec_exclusions, (dot_vec==dot(s1,s2)) )...
    +(R<=0.5).*dot(prob_vec_inclusions, (dot_vec==dot(s1,s2)) );
end

