%% TMz
close all

lambda = 680/10^9;  % wavelength
beta = 2*pi/lambda; % wave vector
a = 0.05*lambda;     % radius of the cylinder
x = -5:0.1:5;
x = x*lambda;
y = x;
[X,Y] = meshgrid(x,y);
[theta,rho] = cart2pol(X,Y);    % Convert cartesian coordinates to polar coordinates
Ei = exp(-j*beta.*rho.*cos(theta)); % incident wave
Es = ScaFeil(beta,a,20,theta,rho);  

figure
I = (abs(Es+Ei)).^2 - (abs(Ei)).^2;
imagesc(I)  % interference field
axis off
axis square
colormap(jet)
colorbar
figure
imagesc((abs(Es)).^2);    % scatter field
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

function Es = ScaFeil(beta,a,n,theta,rho)
    Es = zeros(size(theta));
    Jn = @(nu,x)besselj(nu,x);  % bessel function of first kind
    Hn = @(nu,x)besselh(nu,2,x);    % bessel function of third kind
    x = beta*a;
    for ii = -floor(n/2):ceil(n/2)
        Es = Es + i^(-ii)*Jn(ii,x)/Hn(ii,x)*Hn(ii,beta.*rho).*exp(i*ii*theta);
    end
    ind = find(rho<=a);
    Es(ind) = 0;    % The internal electric field strength of a perfect conductor is 0
end














