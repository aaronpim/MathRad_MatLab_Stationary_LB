function [output] = matrix_index(value,N)
if sign(value)==0
    output = 1:N;
elseif sign(value)==1
    output = 2:N;
else
    output = 1:(N-1);
end

