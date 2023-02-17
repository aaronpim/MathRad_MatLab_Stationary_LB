function [Sctr] = Sigma_s(x,s1,s2,Fokker_Plank)
[X,Y,Z] = meshgrid(x,x,x);
dot_vec = [-1,-2/sqrt(6),-1/sqrt(2),-1/sqrt(3),-1/2,-1/3,0,1/3,1/2,1/sqrt(3),1/sqrt(2),2/sqrt(6),1];
if Fokker_Plank
    ds = norm([1,0,0]-[1,1,0]./sqrt(2));
    if norm(s1-s2)==0
        Sctr = -(4/ds.^2);
    elseif norm(s1-s2)==ds
        Sctr = (1/ds.^2);
    else
        Sctr = 0;
    end
else
    R = sqrt(X.^2+Y.^2+Z.^2);
    prob_vec_inclusions = ones(length(dot_vec),1)/length(dot_vec);
    prob_vec_exclusions = [zeros(1,11),0.1,0.9];
    Sctr = (R>0.5).*dot(prob_vec_exclusions, (dot_vec==dot(s1,s2)) )...
    +(R<=0.5).*dot(prob_vec_inclusions, (dot_vec==dot(s1,s2)) );
end
end

