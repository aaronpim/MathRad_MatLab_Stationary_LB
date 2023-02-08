function [g] = inflow_boundary_data(f,x,S_vec)
%% This function takes an NxNxNxS array, denoted f and imposes Dirichlet 
% boundary conditions on the "inflow".
    [X1,X2]=meshgrid(x,x);
    g = f;
    N = length(x);    
    for i = 1:length(S_vec)
        %If S(1) > 0 then BC on x= -1;
        %If S(1) < 0 then BC on x= +1;
        if sign(S_vec(i,1))>0
            %g(1,:,:,i) = 0.01*exp(-(X1.^2+X2.^2));
            g(1,:,:,i) = 0.1*(sqrt(X1.^2+X2.^2)<0.2).*(sign(S_vec(i,2))==0).*(sign(S_vec(i,3))==0)...
                +0.1*(sqrt(X1.^2+(X2-0.7).^2)<0.2).*(sign(S_vec(i,2))==-1).*(sign(S_vec(i,3))==0)...
                +0.1*(sqrt(X1.^2+(X2+0.7).^2)<0.2).*(sign(S_vec(i,2))==1).*(sign(S_vec(i,3))==0);
        elseif sign(S_vec(i,1))<0
            g(N,:,:,i) = 0;
        end
        %If S(2) > 0 then BC on y= -1;
        %If S(2) < 0 then BC on y= +1;
        if sign(S_vec(i,2))>0
            g(:,1,:,i) = 0;
        elseif sign(S_vec(i,2))<0
            g(:,N,:,i) = 0;
        end
        %If S(3) > 0 then BC on z= -1;
        %If S(3) < 0 then BC on z= +1;
        if sign(S_vec(i,3))>0
            g(:,1,:,i) = 0;
        elseif sign(S_vec(i,3))<0
            g(:,N,:,i) = 0;
        end
        
    end
end

