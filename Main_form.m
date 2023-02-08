%% Parameters of the physical space.
% We consider a cube [-1,1]^3 and discretise it via N. We then conduct
% n_iter iterations.
N = 61;
n_iter = 200;
x = linspace(-1,1,N);
dx = x(2)-x(1);
[Xmat,Ymat] = meshgrid(x,x);
[X,Y,Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);
s_vec = [1,0,0;-1, 0, 0;...
         0,1,0; 0,-1, 0;...
         0,0,1; 0, 0,-1];
%% Cross-Sections
%Sigma_t and Sigma_s are arrays that define the total and scattering 
%cross-sections, which are assumed to be spatially dependent. 
Sigma_t = 0.50*(R<0.5) + 0.90*(R>=0.5);

Sigma_s = @(i,j) (0.95*(R>=0.5) + 0.80*(R<0.5)).*(dot(s_vec(i,:),s_vec(j,:))==1)...
    +(0.04*(R>=0.5) + 0.15*(R<0.5)).*(dot(s_vec(i,:),s_vec(j,:))==0)...
    +(0.01*(R>=0.5) + 0.05*(R<0.5)).*(dot(s_vec(i,:),s_vec(j,:))==-1);

%% Inflow boundary conditions
% The functions g1,...,g6 define the inflow Dirichlet boundary conditions
% on each of the faces.
g1 = 0.01*(sqrt(Xmat.^2+Ymat.^2)<0.6); %0.01*ones(size(Xmat)); %0.01*exp(-sqrt(Xmat.^2+Ymat.^2));
g2 = 0; 
g3 = 0;
g4 = 0;
g5 = 0;
g6 = 0;
%% Array definition
% The array f denotes the solution to the Boltzmann transport equation.
f = zeros(N,N,N,length(s_vec));
f(1,:,:,1) = g1;
f(N,:,:,2) = g2;
f(:,1,:,3) = g3;
f(:,N,:,4) = g4;
f(:,:,1,5) = g5;
f(:,:,N,6) = g6;
f_temp = zeros(N,N,N,length(s_vec));
%% Inflow boundary conditions
for n = 1:n_iter
    for i = 1:length(s_vec)
        sctr = 0;
        for j = 1:length(s_vec)
            sctr = sctr + Sigma_s(i,j).*f(:,:,:,j);
        end
        col = sctr - Sigma_t.*f(:,:,:,i);
        if i == 1
            f_temp(2:N,:,:,1) = f(1:(N-1),:,:,1)+dx*col(1:(N-1),:,:);
            f_temp(1,:,:,1) = g1;
            
        elseif i == 2
            f_temp(1:(N-1),:,:,2) = f(2:N,:,:,2)+dx*col(2:N,:,:);
            f_temp(N,:,:,2) = g2;
            
        elseif i == 3
            f_temp(:,2:N,:,3) = f(:,1:(N-1),:,3)+dx*col(:,1:(N-1),:);
            f_temp(:,1,:,3) = g3;
            
        elseif i == 4
            f_temp(:,1:(N-1),:,4) = f(:,2:N,:,4)+dx*col(:,2:N,:);
            f_temp(:,N,:,4) = g4;
            
        elseif i == 5
            f_temp(:,:,2:N,5) = f(:,:,1:(N-1),5)+dx*col(:,:,1:(N-1));
            f_temp(:,:,1,5) = g5;
            
        else
            f_temp(:,:,1:(N-1),6) = f(:,:,2:N,6)+dx*col(:,:,2:N);
            f_temp(:,:,N,6) = g6;
        end
    end
    f = f_temp;
end
u = sum(f,4);
pcolor(Xmat,Ymat,u(:,:,round(end/2)));
figure;
contourslice(X,Y,Z,u,[-0.75,-0.5,-0.25,0,0.25,0.5,0.75],[],[-0.75,-0.5,-0.25,0,0.25,0.5,0.75])