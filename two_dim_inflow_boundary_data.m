function [g] = two_dim_inflow_boundary_data(f,x,S_vec)
    g = f;
    N = length(x);
    pencil_beam = @(intensity,c,r) intensity.*(abs(x-c)<r);
    index = 1:length(S_vec);
    index_x_equals_minus_1 = index(S_vec(:,1)>0);
    index_x_equals_plus_1 = index(S_vec(:,1)<0);
    index_y_equals_minus_1 = index(S_vec(:,2)>0);
    index_y_equals_plus_1 = index(S_vec(:,2)<0);
%% The boundary x = -1
    for i = index_x_equals_minus_1
        if S_vec(i,2)==0
            g(1,:,i) = pencil_beam(1,0,0.3);
        else
            g(1,:,i) = 0;
            if abs(S_vec(i,1))==2
                g(2,:,i) = 0;
            end
        end
    end
%% The boundary x = 1
    for i = index_x_equals_plus_1
        if S_vec(i,2)==2
            g(N,:,i) = pencil_beam(1,-0.7,0.3);
            g(N-1,:,i) = pencil_beam(1,-0.7,0.3);
        elseif S_vec(i,2)==-1
            if S_vec(i,1)== -1
                g(N,:,i) = pencil_beam(1,0.7,0.2);
            else
                g(N,:,i) = 0;
            end
        else
            g(N,:,i) = 0;
            if abs(S_vec(i,1))==2
                g(N-1,:,i) = 0;
            end
        end
    end
%% The boundary y = -1 
    for i = index_y_equals_minus_1
        g(:,1,i) = 0;
        if abs(S_vec(i,1))==2
            g(:,2,i) = 0;
        end
    end
%% The boundary y = 1 
    for i = index_y_equals_plus_1
        g(:,N,i) = 0;
        if abs(S_vec(i,1))==2
            g(:,N-1,i) = 0;
        end
    end
end

