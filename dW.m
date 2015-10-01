function [ dw ] = dW(q, h)
%dW : Computes derivative of Cubic Spline
%   Detailed explanation goes here

% Ch Ensures normalization
Ch = 5/(14*pi*h*h);
dw = 0;
if q < 1
    dw = Ch * ( -3*((2 - q)^2) + 12*((1-q)^2) ); 
elseif(q < 2)
    dw = Ch * -3 * ((2 - q)^2);
end

end

