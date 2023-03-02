function [Sctr] = Sigma_s(x,s1,s2)
% Input: x % An vector of length N which represents the side length of the
% cuboid domain.

% Input: s1 % A 3x1 vector which represents the trajectory of the particle before scattering

% Input: s2 % A 3x1 vector which represents the trajectory of the particle after scattering

%Output: Sctr % An NxNxN tensor which is such that Sctr(i,j,k) represents
%the scattering cross section sigma_s at the point (x(i),x(j),x(k)) for a
%particle with initial trajectory s1 and scattered trajectory s2.

[X,Y,Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);

dot_vec = [-1,-2/sqrt(6),-1/sqrt(2),-1/sqrt(3),-1/2,-1/3,0,1/3,1/2,1/sqrt(3),1/sqrt(2),2/sqrt(6),1];
% We assume that Sigma_s is a function of dot(s1,s2). Therefore for each
% possible value of the dot product we assign a probability of scattering.

m = length(dot_vec);
prob_vec_inclusions = (m+1)-(1:m); 
prob_vec_inclusions = prob_vec_inclusions./sum(prob_vec_inclusions);

prob_vec_exclusions = 1:m; 
prob_vec_exclusions = prob_vec_exclusions./sum(prob_vec_exclusions);

Sctr = (R>0.5).*dot(prob_vec_exclusions, (abs(dot_vec-dot(s1,s2)))<1e-8)...
    +(R<=0.5).*dot(prob_vec_inclusions, (abs(dot_vec-dot(s1,s2)))<1e-8);
end

