function [g] = inflow_boundary_data(f,x,S_vec)
%% This function takes an NxNxNxS array, denoted f and imposes Dirichlet 
% boundary conditions on the "inflow".
    [X1,X2]=meshgrid(x,x);
    pencil_beam = @(intensity,c,r) intensity.*(sqrt((X1-c(1)).^2+(X2-c(2)).^2)<r);
    g = f;
    N = length(x);
% There are 6 boundaries to consider, and we shall only be interested in
% the trajectory vectors which point "into" the domain. On each boundary
% there are 9 possible trajectories.
index = 1:length(S_vec);
% For each boundary of the cube, we calculate which trajectory vectors are
% "inward pointing"
index_x_equals_minus_1 = index(S_vec(:,1)>0);
index_x_equals_plus_1 = index(S_vec(:,1)<0);
index_y_equals_minus_1 = index(S_vec(:,2)>0);
index_y_equals_plus_1 = index(S_vec(:,2)<0);
index_z_equals_minus_1 = index(S_vec(:,3)>0);
index_z_equals_plus_1 = index(S_vec(:,3)<0);
%% The boundary x = -1
    for i = index_x_equals_minus_1
        % This is an example of a single pencil beam, with radius 0.1,
        % intensity 1, entering in at the point (-1,0,0), and with trajectory
        % (1,0,0); We do not wish for any other beams to enter through this
        % boundary.
        if S_vec(i,2)==0 && S_vec(i,3)==0
            g(1,:,:,i) = pencil_beam(1,[0,0],0.1);
        else
            g(1,:,:,i) = 0;
        end
    end
%% The boundary x = 1
    for i = index_x_equals_plus_1
        % This is an example of a uniform radiance passing through in all
        % possible directions, through the boundary x = 1.
        g(N,:,:,i) = 0.01;
    end
%% The boundary y = -1 
    for i = index_y_equals_minus_1
        % This is an example of a two pencil beama with:
        % radii 0.05 and 0.1,
        % intensities 0.5 and 0.7
        % entry points (0.9,-1,0.9) and (-0.8,-1,0.5)
        % trajectories (-1,1,-1) and (-1,1,0)
        if S_vec(i,1)<0 && S_vec(i,3)<0
            g(:,1,:,i) = pencil_beam(0.5,[0.9,0.9],0.05);
        elseif S_vec(i,1)<0 && S_vec(i,3)==0
            g(:,1,:,i) = pencil_beam(0.7,[-0.8,0.5],0.1);
        else
            g(:,1,:,i) = 0;
        end
    end
%% The boundary y = 1 
    for i = index_y_equals_plus_1
        % If there is no inflow on a boundary then we must explicitly set
        % these values to be zero.
        g(:,N,:,i) = 0;
    end
%% The boundary z = -1 
    for i = index_z_equals_minus_1
        g(:,:,1,i) = 0;
    end
%% The boundary z = 1 
    for i = index_z_equals_plus_1
        g(:,:,N,i) = 0;
    end
end

