function [Sctr] = Sigma_t(x)
% Input: x % A vector of length N which represents the side length of the
% cuboid domain.

%Output: Sctr % An NxN tensor which is such that Sctr(i,j) represents
%the total cross section sigma_t at the point (x(i),x(j)).
    [X,Y] = meshgrid(x,x);
    R = sqrt(X.^2+Y.^2);

    %In this example, we have an inclusion of radius 0.5 which we assume to
    %be denser than the surrounding material.
    Sctr = (R>0.5).*0.6+(R<=0.5).*1.1;
end