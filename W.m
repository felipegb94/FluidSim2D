function [ w ] = W(q, h)
%W : Computes Cubic Spline
%   Detailed explanation goes here

% Ch Ensures normalization
Ch = 5/(14*pi*h*h);
w = 0;
if(1 <= q && q <= 2)
    w = Ch * ((2 - q)^3);
elseif(0 <= q && q < 1)
    w = Ch * ( ((2 - q)^3) - 4*((1-q)^3) );
end

end

