function [g] = two_dim_inflow_boundary_data(f,x,S_vec)
    g = f;
    N = length(x);
    pencil_beam = @(intensity,c,r) intensity.*(abs(x-c)<r);
    index = 1:length(S_vec);
    index_x_equals_minus_1 = index(S_vec(:,1)>0);
    index_x_equals_plus_1 = index(S_vec(:,1)<0);
    index_y_equals_minus_1 = index(S_vec(:,2)>0);
    index_y_equals_plus_1 = index(S_vec(:,2)<0);
%% The boundary y = -1
    for i = index_x_equals_minus_1
        for j=1:S_vec(i,1)
            % This corresponds to a single large pencil beam
            g(j,:,i) = pencil_beam(1,0,0.3).*(S_vec(i,1)==1).*(S_vec(i,2)==0);
        end
    end
%% The boundary y = 1
    for i = index_x_equals_plus_1
        for j=S_vec(i,1):-1
            % This corresponds to no inflow.
            g(N+1+j,:,i) = 0;
        end
    end
    for i = index_x_equals_minus_1
        for j=1:S_vec(i,1)
            % This corresponds to a reflective boundary condition.
            g(N+1-j,:,(S_vec(:,1)==-S_vec(i,1))&(S_vec(:,2)==S_vec(i,2))) = f(N+1-j,:,i);
        end
    end
%% The boundary x = -1 
    for i = index_y_equals_minus_1
        for j=1:S_vec(i,2)
            % This corresponds to no inflow.
            g(:,j,i) = 0;
        end
    end
    for i = index_y_equals_plus_1
        for j=S_vec(i,2):-1
            % This corresponds to a reflective boundary condition.
            g(:,-j,(S_vec(:,1)==S_vec(i,1))&(S_vec(:,2)==-S_vec(i,2))) = f(:,-j,i);
        end
    end
    % This corresponds to a reflective condition.
%% The boundary x = 1
    for i = index_y_equals_plus_1
        for j=S_vec(i,2):-1
            % This corresponds to two pencil beams firing at different
            % angles.
            g(:,N+1+j,i) =   pencil_beam(1,0.0,0.2).*(S_vec(i,1)== 1).*(S_vec(i,2)==-1)...
                            +pencil_beam(1,0.7,0.2).*(S_vec(i,1)==-1).*(S_vec(i,2)==-2);
        end
    end
end

