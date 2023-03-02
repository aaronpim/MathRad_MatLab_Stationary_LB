%% Parameters of the physical space.
% In this section we define the parameters of the physical space by
% considering a cuboid domain [-1,1]^3 and then subdivide it into (N-1)^3
% subcubes. We then iterate the foward Euler process n_iter times.

N = 101;         %Spatial discretisation
n_iter = 600;   %Maximum number of iterations
x = linspace(-1,1,N);           %1-D domain
                 %Spatial step

%% List of trajectory vectors.

s_vec = [-1, 0; 1, 0; 0,-1; 0, 1;...
         -1, 1; 1,-1;-1,-1; 1, 1;...
         -2, 1; 2,-1;-2,-1; 2, 1;...
         -1, 2; 1,-2;-1,-2; 1, 2];
dx = (x(2)-x(1)).*sqrt(sum(s_vec.^2,2));

%% Array definition
% The array f denotes the solution to the Boltzmann transport equation.
f = zeros(N,N,length(s_vec));
f_temp = zeros(N,N,length(s_vec));
f = two_dim_inflow_boundary_data(f,x,s_vec);
total_dose = 0;
count = 0;
%% Finite difference method
while abs(sum(f,'all')-total_dose)>1.0e-6 && count < n_iter
    count = count +1;
    total_dose = sum(f,'all');
    for i = 1:length(s_vec)
        sctr = 0;        
        for j = 1:length(s_vec)
            sctr = sctr + Sigma_s(x,s_vec(i,:),s_vec(j,:)).*f(:,:,j)...
                -Sigma_s(x,s_vec(i,:),s_vec(j,:)).*f(:,:,i);
            %All particles flowing into the point x, minus the particles
            %flowing out of the point x. This is the collision step
        end
        col = sctr - Sigma_t(x).*f(:,:,i); 
        % This step calculates the collision term.
        a1 = matrix_index(s_vec(i,1),N); b1 = matrix_index(-s_vec(i,1),N);
        a2 = matrix_index(s_vec(i,2),N); b2 = matrix_index(-s_vec(i,2),N);
        % This step calculates the indices needed for the finite difference
        % forward Euler method. f(x+ dx*s) = f(x) + dx*col(x)
        f_temp(a1,a2,i) = f(b1,b2,i)+dx(i)*sctr(b1,b2);
    end
    % The temporary vector is to ensure that nothing changes as we iterate
    % throught the trajectory vectors.
    f = f_temp;
    % We impose the boundary conditions again if they have been changed by
    % the finite difference scheme.
    f = two_dim_inflow_boundary_data(f,x,s_vec);
end
u = sum(f,3);
[X,Y] = meshgrid(x,x);
%% Plots
pcolor(X,Y,u);
hold on
plot(0.5*cos(linspace(0,2*pi,101))+((x(2)-x(1))/2),0.5*sin(linspace(0,2*pi,101))+((x(2)-x(1))/2),'r-')
hold off
axis equal
set(gca,'ColorScale','log')