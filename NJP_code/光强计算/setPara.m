function [ para ] = setPara( varargin )
%SSM Summary of this function goes here
% 平板SPR模型，参考论文Principles of nanoparticle imaging using surface plasmons doi:10.1088/1367-2630/17/1/013041
% SetPara: 定义模型的各种变量，可改变的是定义量，生成量根据定义量改变，
% 因此每次改变定义量后，需要刷新生成量，为了效率可以提取一部分刷新，但是有风险。
% 
% 定义量
para.A=1;%系数，与场的入射强度有关
para.dm=47e-9;%金属层厚度
para.e0=8.854187817E-12;%真空介电常数F/m
para.mu0=4*pi*1E-7;%真空磁导率H/m
para.n1=1.514;%玻璃折射率
para.n2=0.1355+1i*3.88;%金膜在680nm下的折射率%0.14737+1i*4.7414;%金膜在780nm下的折射率
para.n3=1.333;%环境介质折射率
para.thetaInc=deg2rad(69.336);%电磁波入射角69.336
para.phiInc=deg2rad(0);%入射方位角
para.wl=670e-9;%电磁波波长[m]
para.rp=50e-9;%颗粒半径[m]
para.np=0.135+3.88*1i;%颗粒折射率
para.siz=[501,501];%模拟区域大小
para.deltaDistance=1e-8;%每像素点间间距[m]
para.center=[251,251];%颗粒中心
para.nLimit=5;%求和上限
para.deltah=200e-9;%颗粒距金膜高度
% 定义量修改

while length(varargin)>=2
    prop =varargin{1};
    val=varargin{2};
    varargin=varargin(3:end);
    eval(['para.' prop '=val;']);     
end

% 生成量
para.h=para.rp+para.deltah;%颗粒中心坐标
para.z=-para.dm;%视野距金膜下表面距离
para.k0=2*pi/para.wl;%波矢大小
para.e3=para.n3^2;%环境介电常数
para.e1=para.n1^2;
para.omega=2*pi*3e8/para.wl;%圆频率
% [X,Y]=meshgrid(1:para.siz(1),1:para.siz(2));
[X,Y]=meshgrid(1:para.siz(2),1:para.siz(1));
para.r=sqrt((para.deltaDistance*sqrt((X-para.center(2)).^2+(Y-para.center(1)).^2)).^2+para.h.^2);%图上点的球坐标r的坐标
para.theta=acos(para.h./para.r);%天顶角[0,pi/2)%asin((X-para.center(2))./para.r*para.deltaDistance)-pi/2;%图上点的球坐标theta的坐标[0,2pi)
%para.phi=0;%图上点的球坐标theta的坐标
para.phi=zeros(para.siz(1),para.siz(2));
for iy=1:para.siz(1)
    for ix=1:para.siz(2)
        dx=ix-para.center(2);
        dy=iy-para.center(1);
        if dx>=0 && dy>=0
            para.phi(iy,ix)=atan(dy/dx);
        end
        if dx<0 && dy>=0
            para.phi(iy,ix)=pi-atan(-dy/dx);
        end
        if dx<0 && dy<0
            para.phi(iy,ix)=pi+atan(dy/dx);
        end
        if dx>=0 && dy<0
            para.phi(iy,ix)=2*pi-atan(-dy/dx);
        end%方位角phi=0-2pi
    end
end

para.kx=para.n1*abs(para.k0)*sin(para.thetaInc)*cos(para.phiInc);
para.ky=para.n1*abs(para.k0)*sin(para.thetaInc)*sin(para.phiInc);

% para.k3z=-para.n3*abs(para.k0)*sqrt(1-(para.n1/para.n3*sin(para.thetaInc)));
% para.k2z=para.k3z.*para.n2./para.n3;
% para.k1z=para.k3z.*para.n1./para.n3;
% % 以上三个波矢z分量原作此，今依论文更正
para.k3z=-para.n3*abs(para.k0)*sqrt(1-(para.n1/para.n3*sin(para.thetaInc))^2);
para.k2z=para.n2*abs(para.k0)*sqrt(1-(para.n1/para.n2*sin(para.thetaInc))^2);
para.k1z=para.n1*abs(para.k0)*cos(para.thetaInc);

para.kxy=sqrt(para.kx^2+para.ky^2);
para.N=para.np/para.n3;
para.rho=para.kxy*para.rp;

para.B=para.A.*(exp((abs(para.k3z)-abs(para.k1z)).*para.dm./2)./(2.*abs(para.k2z)./(para.n2.*para.n2))).*(exp(-abs(para.k2z).*para.dm).*((abs(para.k2z)./(para.n2.*para.n2))-(abs(para.k1z)./(para.n1.*para.n1)))+exp(abs(para.k2z).*para.dm).*((abs(para.k1z)./(para.n1.*para.n1))+(abs(para.k2z)./(para.n2.*para.n2))));%系数B
para.Ii=abs(para.B).*abs(para.B).*para.e3.*para.e3.*pi.*sinh(2.*abs(para.k3z).*para.rp).*exp(abs(para.k3z).*para.h)./(abs(para.k3z).*para.rp);%SPP波入射到纳米颗粒上的强度

para.E0x=para.B.*exp(abs(para.k3z).*para.z)./(para.omega.*para.e3.*para.e0).*1i.*para.k3z.*...
    exp(1i*para.deltaDistance*(para.kx*(X-para.center(2))+para.ky*(Y-para.center(1)))); % 即正文式1
for n=0:para.nLimit
    [para.ANR(n+1),para.ANTHETA(n+1)]=an(n,para.N,para.rho);
end

% para.E0r=ones(para.siz);
% 此句极怪，下面依文章更订
para.E0r = -i*abs(para.k3z).*sin(para.theta).*cos(para.phi) ...
    - i*para.ky*abs(para.k3z)/para.kx.*sin(para.theta).*sin(para.phi)...
    + (para.kx + (para.ky)^2/para.kx)*cos(theta);
para.E0theta = -i*abs(para.k3z).*cos(para.theta).*cos(para.phi) ...
    - i*para.ky*abs(para.k3z)/para.kx.*cos(para.theta).*sin(para.phi)...
    + (para.kx + (para.ky)^2/para.kx)*sin(theta);
para.E0phi = i*abs(para.k3z)*sin(para.phi) - i*para.ky*abs(para.k3z)/para.kx*cos(para.phi);

end

