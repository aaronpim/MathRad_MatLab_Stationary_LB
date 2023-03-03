function [Sctr] = Sigma_s(x,s1,s2,s_angle)
%Input: x % An vector of length N which represents the side length of the
% square domain.

% Input: s1 % A 2x1 vector which represents the trajectory of the particle before scattering

% Input: s2 % A 2x1 vector which represents the trajectory of the particle after scattering

% Input: s_angle % The set of all possible angles between two vectors in
% the discrete set of trajectory vectors. This vector is assumed to be
% ordered lowest angle to highest e.g. -pi to pi.

%Output: Sctr % An NxN tensor which is such that Sctr(i,j) represents
%the scattering cross section sigma_s at the point (x(i),x(j)) for a
%particle with initial trajectory s1 and scattered trajectory s2.

    [X,Y] = meshgrid(x,x);
    R = sqrt(X.^2+Y.^2);
    N = length(s_angle);
    
    %This represents the scattering distribution in the inclusion.
    prob_vec_inclusions = ones(N,1); 
    prob_vec_inclusions = prob_vec_inclusions./sum(prob_vec_inclusions);
    %The above represents a uniform distribution in scattering, each angle is
    %equally likely

    %This represents the scattering distribution in the bulk domain.
    prob_vec_exclusions = [zeros(N-1,1);1]; 
    prob_vec_exclusions = prob_vec_exclusions./sum(prob_vec_exclusions);
    %The above represents a dirac delta in terms of trajectory, there is
    %scattering in only one direction, which is directly forwards. i.e. no
    %scattering.

    %We define the Bulk domain to be the square outside of a sphere of
    %radius 0.5 and the Inclusion being inside the sphere of radius 0.5.
    [~,I] = min(s_angle-subspace(s1',s2'));
    %We use min(theta-arccos(s.s')) rather than theta == arccos(s.s') to 
    %account for machine error.
    Sctr = (R>0.5).*prob_vec_exclusions(I)+(R<=0.5).*prob_vec_inclusions(I);
end

