function [f,not_holes] = domain_holes(f,x)
%% This function defines holes in the domain
[X,Y]=meshgrid(x,x);

% This generates a rectangular hole with centre (0.7,0.4), width 0.2 and
% height 0.1
rct_cntr = [0.7,0.4]; rct_wdth = 0.2; rct_hth = 0.1;
square_hole = (2*abs(X-rct_cntr(1))<=rct_wdth)&(2*abs(Y-rct_cntr(2))<=rct_hth);

% This generates a circular hole with centre (0.1,-0.5) and radius 0.25
crc_cntr = [0.1,-0.5]; crc_rad = 0.25;
circle_hole = (sqrt((X-crc_cntr(1)).^2+(Y-crc_cntr(2)).^2)<crc_rad);

%We then combine the holes and define the part of the domain that is not a
%hole. The variable not_holes(i,j) is such that it is 1 if the point
%(x(i),x(j)) is not a hole and 0 if it is a hole.
not_holes = ~ (square_hole | circle_hole);

%For each possible trajectory, if you are not a hole you retain your value,
%otherwise you are set to zero.
for i = 1:size(f,3)
    f(:,:,i)=f(:,:,i).*not_holes;
end
end