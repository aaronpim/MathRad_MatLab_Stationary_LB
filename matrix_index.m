function [output] = matrix_index(value,N)
if value == 0
    output = 1:N;
elseif value ==1
    output = 2:N;
elseif value ==2
    output = 3:N;
elseif value ==-1
    output = 1:(N-1);
elseif value ==-2
    output = 1:(N-2);
end