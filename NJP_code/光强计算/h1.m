function [ h ] = h1( n,x )
%H1 Summary of this function goes here
%   球化的第三类bessel函数

h=sqrt(pi/2./x).*besselh(n+0.5,1,x);

end

