function [ dw ] = dW(q, h)
%dW : Computes derivative of Cubic Spline
%   Detailed explanation goes here

% Ch Ensures normalization
Ch = 15/(14*pi*h*h);
dw = 0;
if(1 <= q && q <= 2)
    dw = Ch * -3 * ((2 - q)^2);
elseif(0 <= q && q < 1)
    dw = Ch * ( -3*((2 - q)^2) + 12*((1-q)^2) );
end

end

