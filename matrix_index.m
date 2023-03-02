function [output] = matrix_index(value,N)
%This function controls the indexing of the tensor for the Euler method
% dot(s,grad(f)) ~= (f(x+dx*s)-f(x))/dx

% f(x+dx*s) is represented by
% f(matrix_index(s(1)),matrix_index(s(2)),matrix_index(s(3)))

%Similiarly f(x) is reperesented by 
% f(matrix_index(-s(1)),matrix_index(-s(2)),matrix_index(-s(3)))
if sign(value)==0
    output = 1:N;
elseif sign(value)==1
    output = 2:N;
else
    output = 1:(N-1);
end

