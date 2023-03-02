function [Sctr] = Sigma_t(x)
% Input: x % A vector of length N which represents the side length of the
% cuboid domain.

%Output: Sctr % An NxNxN tensor which is such that Sctr(i,j,k) represents
%the total cross section sigma_t at the point (x(i),x(j),x(k)).

[X,Y,Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);
Sctr = 0.05*(R>0.5)+0.5*(R<=0.5);
end

