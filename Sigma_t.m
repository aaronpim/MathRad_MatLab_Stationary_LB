function [Sctr] = Sigma_t(x)
[X,Y,Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);

Sctr = 1.0*(R>0.5)+1.0*(R<=0.5);
end

