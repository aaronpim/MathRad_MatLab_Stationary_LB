function [Sctr] = Sigma_t(x)
    [X,Y] = meshgrid(x,x);
    R = sqrt(X.^2+Y.^2);
    Sctr = (R>0.5).*0.6+(R<=0.5).*1.1;
end