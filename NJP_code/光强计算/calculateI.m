%计算光强I
%2020.08.06

%%
[para]=setPara( );
[E_image,E_rad_s,E_spp_s]=E_s(para);
% I=abs(para.E0x+E_rad_s+E_spp_s).^2-abs(para.E0x).^2;
IM=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
figure;imagesc(IM);

%%
[para1]=setPara('np',1.373);
[E_image1,E_rad_s1,E_spp_s1]=E_s(para1);
% I1=abs(para1.E0x+E_rad_s1+E_spp_s1).^2-abs(para1.E0x).^2;
Ii1=abs(para1.E0x+E_image1).^2-abs(para1.E0x).^2;

% [para2]=setPara('np',1.5);
% [E_image2,E_rad_s2,E_spp_s2]=E_s(para2);
% I2=abs(para2.E0x+E_rad_s2+E_spp_s2).^2;
% Ii2=abs(para2.E0x+E_image2).^2;

[para2]=setPara('np',1.37);
[E_image2,E_rad_s2,E_spp_s2]=E_s(para2);
% I2=abs(para2.E0x+E_rad_s2+E_spp_s2).^2-abs(para2.E0x).^2;
Ii2=abs(para2.E0x+E_image2).^2-abs(para2.E0x).^2;

% [para3]=setPara('np',1.8);
% [E_image3,E_rad_s3,E_spp_s3]=E_s(para3);
% I3=abs(para3.E0x+E_rad_s3+E_spp_s3).^2;
% Ii3=abs(para3.E0x+E_image3).^2;

%%
np=1.45:0.0001:1.48;
for i=1:length(np)
    para=setPara('np',np(i));
    [E_image,E_rad_s,E_spp_s]=E_s(para);
    Ii=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
    Itotal(i,:,:)=Ii;
    imagetotal(i,:,:)=E_image;
    radtotal(i,:,:)=E_rad_s;
    spptotal(i,:,:)=E_spp_s;
    paratotal(i)=para;
    disp(['np=',num2str(np(i)),',timenow:',datestr(now)]);
end
showSlide(Itotal);
%%
for i=1:(length(Itotal(:,1,1)))
    per(i,:,:)=(Itotal(i,:,:)-Itotal(1,:,:))./Itotal(1,:,:);
end
showSlide(per,'elsecommand','caxis([0,0.1])');
%%
for i=1:(length(Itotal(:,1,1)))
    dif(i,:,:)=(Itotal(i,:,:)-Itotal(1,:,:));
end
showSlide(dif,'elsecommand','caxis([0,1e-8])');
%%
roi1=squeeze(Itotal(1,230:270,385:495));
num1=sum(sum(roi1));
for i=1:(length(Itotal(:,1,1)))
    roi2=squeeze(Itotal(i,230:270,385:495));
    num2=sum(sum(roi2));
    result(i)=(num2-num1)/num1;
end
%%
for i=1:31
    roi(i,:,:)=squeeze(Itotal(i,230:270,385:495));
    num(i)=sum(sum(squeeze(roi(i,:,:))))./(41*111);
    result(i)=num(i)/(2e-7);
end
%%
load('theta.mat');
clear lambda;

load('nkans.mat');
lambda=(nkans(:,1))';ngold=(nkans(:,6))';kgold=(nkans(:,7))';

wl=680e-9:5e-9:785e-9;
ilambda=(680-299):5:(785-299);
itheta=(680-399):5:(785-399);

for i=1:length(wl)
    n=ngold(ilambda(i))+1i*kgold(ilambda(i));
    anstheta=deg2rad(theta(itheta(i)))+0.5*pi/180;
    [para]=setPara('n2',n,'np',n,'thetaInc',anstheta,'wl',wl(i));
    [E_image,E_rad_s,E_spp_s]=E_s(para);
    I=abs(para.E0x+E_rad_s+E_spp_s).^2;
    Ii=abs(para.E0x+E_image).^2;
    filename=[num2str(wl(i)),'.mat'];
    save(filename,'E_image','E_rad_s','E_spp_s','para','I','Ii');
end
%%


[para1]=setPara('np',1.46);
[E_image1,E_rad_s1,E_spp_s1]=E_s(para1);
Ii1=abs(para1.E0x+E_image1).^2-abs(para1.E0x).^2;

[para2]=setPara('np',1.463);
[E_image2,E_rad_s2,E_spp_s2]=E_s(para2);
Ii2=abs(para2.E0x+E_image2).^2-abs(para2.E0x).^2;


% [para3]=setPara('np',1.8);
% [E_image3,E_rad_s3,E_spp_s3]=E_s(para3);
% Ii3=abs(para3.E0x+E_image3).^2-abs(para3.E0x).^2;

per21=(Ii2-Ii1)./Ii1;
% per31=(Ii3-Ii1)./Ii1;

%%
figure;imagesc(Ii1);
figure;imagesc(Ii2);
% figure;imagesc(Ii3);
figure;imagesc(per21);caxis([0,0.1]);
% figure;imagesc(per31);caxis([0,1]);
%%
% load('1028.mat');
per21=(Ii2-Ii1)./Ii1;
% imagesc(per21);caxis([0,0.1]);
roi1=Ii1(230:270,385:495);
roi2=Ii2(230:270,385:495);
num1=sum(sum(roi1));
num2=sum(sum(roi2));
result=(num2-num1)/num1;
%%
Ii=zeros(40,361,361);
for i=1:40
    rp=i*5e-9;
    h=rp+5e-9;
    para=setPara('rp',rp,'h',h);
    [E_image,E_rad_s,E_spp_s]=E_s(para);
    Ii(i,:,:)=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
end
%%
w=680:1:790;
wl=w.*1e-9;

load('nkans.mat');
ngold=(nkans(:,6))';kgold=(nkans(:,7))';
n2=ngold+1i.*kgold;
cn2=299;%lambda从300开始

load('theta.mat');
thetaInc=deg2rad(theta(:,2));
cthetaInc=399;%从400开始

% load('lambda_nk_ag.mat');
% np=lambda_nk_ag(:,2)+1i.*lambda_nk_ag(:,3);
% cnp=299;
np=n2;cnp=cn2;%Au

for rp=30e-9
    h=rp+1e-9;   
    for i=1:length(w)
        [para]=setPara('wl',wl(i),'n2',n2(w(i)-cn2),'thetaInc',thetaInc(w(i)-cthetaInc),'np',np(w(i)-cnp),'rp',rp,'h',h);
        [E_image,E_rad_s,E_spp_s]=E_s(para);
        Ii=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
        filename=['Au_rp=',num2str(rp),'_w=',num2str(w(i)),'.mat'];
        save(filename,'E_image','E_rad_s','E_spp_s','para','Ii');
    end
end
%% 20210305溶液折射率虚部的变化
n3i=0.871:0.001:0.871;
for i=1:length(n3i)
    n3=1.333+1i*n3i(i);
    para=setPara('n3',n3);
    [E_image,E_rad_s,E_spp_s]=E_s(para);
    filename=['E:\MATLAB\202008_光强计算_来源IOP\20210309\n3i=',num2str(n3i(i)),'.mat'];
    save(filename,'E_image','E_rad_s','E_spp_s','para');
    disp(['n3i=',num2str(n3i(i)),',timenow:',datestr(now)]);
    clear E_image E_rad_s E_spp_s para
end
%% 20210307 one data
n3i=0:0.001:1;
for i=1:length(n3i)
    load(['E:\MATLAB\202008_光强计算_来源IOP\20210309\n3i=',num2str(n3i(i)),'.mat']);
    Ii=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
    totalI(i,:,:)=Ii;
end
%% 20210412
for i=5:10
    rp=i*10e-9;
    h=rp+5e-9;
    para=setPara('rp',rp,'h',h);
    [E_image,E_rad_s,E_spp_s]=E_s(para);
    Ii(i,:,:)=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
end
%%
for i=1:15
    np=1+0.1*i;
    para=setPara('np',np);
    [E_image,E_rad_s,E_spp_s]=E_s(para);
    Ii(i,:,:)=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
end
%% 20211201
rplist=50e-9:10e-9:100e-9;
nplist=0.135+3.88*1i;
auIm=zeros(length(rplist),501,501);
for ir=1:length(rplist)
    in=1;
    [para]=setPara('rp',rplist(ir),'np',nplist(in));
    [E_image,~,~]=E_s(para);
    Ii=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
    auIm(ir,:,:)=Ii;
end
%% 20211202
rplist=50e-9:10e-9:100e-9;
auIm=zeros(length(rplist),501,501);
for ir=1:length(rplist)
    [para]=setPara('rp',rplist(ir));
    [E_image,~,~]=E_s(para);
    Ii=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
    auIm(ir,:,:)=Ii;
end
%% 20220116
thetalist=65:1:75;
Im=zeros(length(thetalist),501,501);
datelist=['Start - timenow:',datestr(now)];
disp(datelist);
TimeText=datelist;
for i=1:length(thetalist)
    [para]=setPara('thetaInc',deg2rad(thetalist(i)));
    [E_image,~,~]=E_s(para);
    Ii=abs(para.E0x+E_image).^2-abs(para.E0x).^2;
    Im(i,:,:)=Ii;
    datelist=['theta: ',num2str(thetalist(i)),'/',num2str(length(thetalist)),' - timenow:',datestr(now)];
    disp(datelist);
    TimeText=[datelist newline TimeText];
end