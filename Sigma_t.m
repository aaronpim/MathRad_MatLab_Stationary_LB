function [Sctr] = Sigma_t(x)
[X,Y,Z] = meshgrid(x,x,x);
R = sqrt(X.^2+Y.^2+Z.^2);

Sctr = 0.95*(R>0.5)+0.9*(R<=0.5);
end

