
theta = 0:0.01*pi:2*pi;
x = cos(theta);
y = cot(theta);
for n = 2:5
    P = legendre(n,x);
    tao_n = y.*P(2,:) - P(3,:);
    figure
    polarplot(theta,abs(tao_n));
end

%% 辅助函数τ和π
close all
theta = 0:0.01*pi:2*pi;
x = cos(theta);
for ii = 1:5
    [pi_n,tao_n] = pi_tao(ii,x);
    figure
    subplot(121)
    polarplot(theta,(tao_n));   % 这里考虑了半径的正负，但是极坐标里面显示不出来
    subplot(122)
    polarplot(theta,(pi_n));    % 参考MATLAB官方文档可以显示负半径
end


function [pi_n,tao_n] = pi_tao(n,x)
    
    pi0 = 0;
    pi1 = 1;
    
    if n == 0
        pi_n = 0;
    elseif n == 1
        pi_n = 1*ones(size(x));
        tao_n = x;
    else
        a = pi1;
        b = pi0;
        for ii = 2:n
            pi_n = (2*ii-1)/(ii-1)*x.*a - ii/(ii-1)*b;
            tao_n = ii*x.*pi_n - (ii+1)*a;
            b = a;
            a = pi_n;
        end
    end
end

%% 计算矢量波函数M和N
lambda = 670/1e9;
N = 1.5;
k = 2*pi/lambda*N;
function [M,N] = VecSphHarm(n,k,az,el,r)
    rho = k*r;
    M = cell(3,1);
    N = M;
    syms jn(n,x)
    jn(n,x) = sqrt(pi/2./x).*besselj(n+0.5,x);
    rhoj(n,x) = x.*jn;
    drhoj = diff(rhoj,x,1);
    jn = matlabFunction(jn);
    rhoj = matlabFunction(rhoj);
    drhoj = matlabFunction(drhoj);

    M{1} = cos(az)*n*(n+1).*sin(el).*pi_n.*jn(rho)/rho;
    M{2} = cos(az).*tao_n.*diff(rhoj(rho))./rho;
    M{3} = -sin(az).*pi_n.*diff(rhoj(rho))./rho;
end

    x = 0:0.1:10;
    x = x*lambda;
    y = x;
    z = x;
    [x,y,z] = meshgrid(x,y,z);
    [az,el,r] = cart2sph(x,y,z);

































