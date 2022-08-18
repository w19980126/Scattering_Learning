%% TMz
close all

epsilonr = 2.55;
lambda = 680/10^9;  % wavelength
k0 = 2*pi/lambda; % wave vector
k1 = sqrt(epsilonr)*k0;
a = 0.05*lambda;     % radius of the cylinder
x = -50:0.2:50;
x = x*lambda;
y = x;
[X,Y] = meshgrid(x,y);
[theta,rho] = cart2pol(X,Y);    % Convert cartesian coordinates to polar coordinates
Ei = exp(-j*k0.*rho.*cos(theta)); % incident wave
[Es ,Et]= ScaFeil(k0,k1,a,21,theta,rho,epsilonr);  

figure
I = (abs(Es+Ei+Et)).^2-(abs(Ei)).^2;
imagesc(I)  % interference field
axis off
axis square
colormap(jet)
colorbar
figure
imagesc((abs(Es)).^2)    % scatter field
axis off
axis square
colormap(jet)
colorbar
F = fftshift(fft2(I));
F(101,101) = 0;
figure
imagesc(abs(F))    % k space
axis off
axis square
colormap(jet)
colorbar

function [Es,Et] = ScaFeil(k0,k1,a,n,theta,rho,epsilonr)
    Es = zeros(size(theta));
    Et = zeros(size(theta));

    mu0 = 4*pi*10^(-7);
    epsilon0 = 8.8542e-12;

    syms Jn(nu,x) Hn(nu,x)
    Jn(nu,x) = besselj(nu,x);  % bessel function of first kind
    dJn = diff(Jn,x);
    Hn(nu,x) = besselh(nu,2,x);    % bessel function of third kind
    dHn = diff(Hn,x);
    x = k0*a;
    Jn = matlabFunction(Jn);
    Hn = matlabFunction(Hn);
    dHn = matlabFunction(dHn);
    dJn = matlabFunction(dJn);

    x = k0*a;
    x1 = k1*a;
    eta = sqrt(mu0/epsilon0);
    eta1 = sqrt(mu0/epsilon0/epsilonr);
    for ii = -floor(n/2):floor(n/2)
        an = i^(-ii)*(dJn(ii,x)*Jn(ii,x1)*eta1-Jn(ii,x)*dJn(ii,x1)*eta)...
            /(Hn(ii,x)*dJn(ii,x1)*eta-dHn(ii,x)*Jn(ii,x1)*eta1);    % scattering coefficient
        bn = i^(-ii)*(Hn(ii,x)*dJn(ii,x)*eta1-Jn(ii,x)*dHn(ii,x)*eta1)...
            /(Hn(ii,x)*dJn(ii,x1)*eta-dHn(ii,x)*Jn(ii,x1)*eta1);    % transmission coefficient
        Es = Es + an*Hn(ii,k0.*rho).*exp(i*ii*theta);
        Et = Et + an*bn*Jn(ii,k1.*rho).*exp(i*ii*theta);
    end
    ind = find(rho<=a);
    Es(ind) = 0;    % The internal electric field strength of a perfect conductor is 0
    ind = find(rho>a);
    Et(ind) = 0;
end
























