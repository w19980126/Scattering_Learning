function [ E_image,E_rad_s,E_spp_s ] = E_s( para )%e.g. [ E_rad_s,E_spp_s ]=E_s(setPara( ));
% 主计算函数，计算模型表面等离子体波E_spp_s，以及脱离表面的电场E_rad_s
    E_rad_s=zeros(size(para.E0r));
    E_spp_s=E_rad_s;
    E_image=E_rad_s;
    parfor N=1:numel(para.r)
        [E_image(N),E_rad_s(N),E_spp_s(N)]=E_s_point( N,para );%可以改进
    end
end


%%
function [ E_image,E_rad_s,E_spp_s ] = E_s_point( N,para )

    tempSR=0;
    tempSTHETA=0;
    for n=0:para.nLimit
        tempR=4*pi*1i^n*para.ANR(n+1)*h1(n,para.kxy*para.r(N))*para.E0r(N);
        tempTheta=4*pi*1i^n*para.ANTHETA(n+1)*h1(n,para.kxy*para.r(N))*para.E0theta(N);
        tempY=0;
        for m=-n:n
            tempY=tempY+Y(n,m,para.theta(N),para.phi(N))*conj(Y(n,m,para.thetaInc,para.phiInc));
        end
        tempSR=tempSR+tempR*tempY;
        tempSTHETA=tempSTHETA+tempTheta*tempY;
    end
    evanescentRp=para.B*exp(abs(para.k3z)*para.rp*cos(para.theta(N)))/(para.omega*para.e3*para.e0);
    E_rad_s=evanescentRp*tempSR;
    
    evanescentImage=para.A.*exp(-abs(para.k1z).*para.dm./2-1i.*abs(para.k1z).*para.z)./(para.omega.*para.e0.*para.e1.*abs(para.np).*sqrt(para.Ii));
    E_image=evanescentImage*tempSR;
    
    evanescentR=para.B*exp(abs(para.k3z)*para.r(N)*cos(para.theta(N)))/(para.omega*para.e3*para.e0);
    E_spp_s=evanescentR*tempSTHETA;

end

