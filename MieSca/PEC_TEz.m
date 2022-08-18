% %% TEz electric field with something wrong
% close all
% 
% lambda = 680/10^9;  % wavelength
% k0 = 2*pi/lambda; % wave vector
% a = 0.5*lambda;     % radius of the cylinder
% 
% x = -5:0.1:5;
% x = x*lambda;
% y = x;
% [X,Y] = meshgrid(x,y);
% [theta,rho] = cart2pol(X,Y);    % Convert cartesian coordinates to polar coordinates
% [EiR,EiP,EsR,EsP] = ScaFeil(k0,a,51,theta,rho);
% 
% figure
% I = (abs(EiR+EsR)).^2;
% imagesc(I)  % interference field
% axis off
% axis square
% figure
% imagesc(abs(EsR+EsP))    % scatter field
% axis off
% axis square
% F = fftshift(fft2(I));
% F(101,101) = 0;
% figure
% imagesc(log(abs(F)))    % k space
% axis off
% axis square
% 
% function [EiR,EiP,EsR,EsP] = ScaFeil(k0,a,n,theta,rho)
%     lambda = 680/10^9;  % wavelength
%     k0 = 2*pi/lambda; % wave vector
%     mu0 = 4*pi*10^(-7);
%     epsilon0 = 8.8542e-12;
%     H0 = 1;
%     omega = k0/sqrt(mu0*eps);
% 
%     eta = sqrt(mu0/epsilon0);
%     EiR = zeros(size(theta));
%     EiP = zeros(size(theta));
%     EsR = zeros(size(theta));
%     EsP = zeros(size(theta));
%     syms Jn(nu,x) Hn(nu,x)
%     Jn(nu,x) = besselj(nu,x);  % bessel function of first kind
%     dJn = diff(Jn,x);
%     Hn(nu,x) = besselh(nu,2,x);    % bessel function of third kind
%     dHn = diff(Hn,x);
%     x = k0*a;
%     Jn = matlabFunction(Jn);
%     Hn = matlabFunction(Hn);
%     dHn = matlabFunction(dHn);
%     dJn = matlabFunction(dJn);
% 
%     for ii = -floor(n/2):ceil(n/2)
%         EiR = EiR + H0/omega/epsilon0./rho*ii*i^(-ii).*Jn(ii,k0.*rho).*exp(i*ii*theta);
%         EiP = EiP + i*eta*H0*i^(-ii).*dJn(ii,k0.*rho).*exp(i*ii*theta);
%         EsR = EsR - H0/omega/epsilon0./rho*ii*i^(-ii)*dJn(ii,x)/dHn(ii,x).*Hn(ii,k0.*rho).*exp(i*ii*theta);
%         EsP = EsP - i*eta*H0*i^(-ii)*dJn(ii,x)/dHn(ii,x).*dHn(ii,k0.*rho).*exp(i*ii*theta);
%     end
% %     ind = find(rho<=a);
% %     Es(ind) = 0;    % The internal electric field strength of a perfect conductor is 0
% end


%% TMz magnetic field
close all
lambda = 680/10^9;  % wavelength
beta = 2*pi/lambda; % wave vector
a = 0.5*lambda;     % radius of the cylinder
x = -10:0.02:10;
x = x*lambda;
y = x;
[X,Y] = meshgrid(x,y);
[theta,rho] = cart2pol(X,Y);    % Convert cartesian coordinates to polar coordinates
Hi = exp(-j*beta.*rho.*cos(theta)); % incident wave
Hs = ScaFeil(beta,a,51,theta,rho);  

figure
I = (abs(Hs+Hi)).^2;
imagesc(I)  % interference field
axis off
axis square
figure
imagesc((abs(Hs)).^2)    % scatter field
axis off
axis square
F = fftshift(fft2(I));
F(101,101) = 0;
figure
imagesc(log(abs(F)))    % k space
axis off
axis square

function Hs = ScaFeil(beta,a,n,theta,rho)
    Hs = zeros(size(theta));
    syms Jn(nu,x) Hn(nu,x)
    Jn(nu,x) = besselj(nu,x);  % bessel function of first kind
    dJn = diff(Jn,x);
    Hn(nu,x) = besselh(nu,2,x);    % bessel function of third kind
    dHn = diff(Hn,x);
    Jn = matlabFunction(Jn);
    Hn = matlabFunction(Hn);
    dHn = matlabFunction(dHn);
    dJn = matlabFunction(dJn);
    x = beta*a;
    
    for ii = -floor(n/2):floor(n/2)
        an = -i^(-ii)*dJn(ii,x)/dHn(ii,x);
        Hs = Hs + an*Hn(ii,beta.*rho).*exp(i*ii*theta);
    end
    ind = find(rho<=a);
    Hs(ind) = 0;    % The internal electric field strength of a perfect conductor is 0
end









