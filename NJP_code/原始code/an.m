function [ ANR,ANTHETA ] = an( n,N,rho )
%ANR Summary of this function goes here
%   Detailed explanation goes here
Nrho=N*rho;
ANTHETA=(Jn(n,Nrho)*dJn(n,rho)-N*Jn(n,rho)*dJn(n,Nrho))/(N*h1(n,rho)*dJn(n,Nrho)-Jn(n,Nrho)*dh1(n,rho));
ANR=-1i*(N-1)/(N+1)*ANTHETA;

end

function [ data ] = dJn( n,x )
%ANR Summary of this function goes here
%   Detailed explanation goes here
data=Jn(n-1,x)-(n+1)/x*Jn(n,x);

end

function [ data ] = dh1( n,x )
%ANR Summary of this function goes here
%   Detailed explanation goes here
data=(h1(n-1,x)-(h1(n,x)+x*h1(n+1,x))/x)/2;

end