%% 通过解调恢复球形散射场
theta = 0;
n = 1.5;
lambda = 680;
kapa = 15;
phi = (0:0.05:2)*pi;
scale_factor = 0.1;
M_size = 401;
for ii = 1:length(phi)
    [Ei,Es,F,I] =  wave_generate(lambda,n,kapa,theta,phi(ii),scale_factor,M_size,90,70.5);
    temp = (I - min(min(I)))/(max(max(I)) - min(min(I)));
    imshow(temp)
    axis off
    axis square
    title([sprintf('Psi = %.2f',phi(ii)/pi),'\pi'],'fontsize',15,'fontweight','bold');
    pause(0.1)
end

mask = zeros(M_size);
K = fftshift(fft2(I.*cos(angle(Es*exp(-i*phi(ii))))));  % 乘以载波后的k空间
imagesc(abs(K));
mask(171:231,171:231) = 1;  % 低通滤波
S = abs(ifft2(ifftshift(K.*mask))); % 重构的散射场
imagesc(S)

%% 计算Psi，原理参考_根据希尔伯特变换恢复散射场PPT_
theta = 0;
n = [1.51,1.33];
lambda = 680;
kapa = 3;
phi = (0:0.05:2)*pi;
scale_factor = 0.1;
M_size = 401;
mask = zeros(M_size);
mask(171:231,171:231) = 1;  % 低通滤波
for ii = 1:length(phi)
    [Ei,Es,F,I] = wave_generate(lambda,n,kapa,theta,phi(ii),scale_factor,M_size);
    K = fftshift(fft2(I.*cos(angle(Es*exp(-i*phi(ii))))));  % 乘以载波后的k空间
    S = abs(ifft2(ifftshift(K.*mask))); % 重构的散射场
    psi(ii) = acos(S(201,201)/abs(Es(201,201)));
    c0(ii) = S(201,201);
end
figure
plot(phi,psi)

%% 希尔伯特变换
plot(real(2*Es(201,:)))
hold on
plot(abs(hilbert(I(201,:))))
plot(abs(Es(201,:))*2)
plot(cos(angle(Es(201,201:end))))
hold on
plot(I(201,:))
plot(abs((fft(I(201,201:end)))))

plot(I(201,201:end)./(abs(hilbert(I(201,201:end)))))
plot(abs(fft(cos(angle(Es(201,201:end))))))

mytest = @(x,xdata) cos(x(1).*xdata+x(2));
r = lsqcurvefit(mytest,[2*pi/5.1538,pi],0:200,I(201,201:end)./abs(hilbert(I(201,201:end))));
figure
plot(cos(angle(Es(201,201:end))))
hold on
plot(1:201,mytest(r,0:200))

plot(1:201,cos(2*pi/5.1538*(0:200)+phi(ii)))

plot(abs(fft(cos(2*pi/5.1538*(1:201)))))

%% 通过希尔伯特变换辅助恢复散射场强度计算Psi值
theta = 0;
n = [1.51,1.33];
lambda = 680;
kapa = 3;   % SPP衰减长度
phi = (0:0.05:2)*pi;    % 预设的Psi值，记为phi
scale_factor = 0.1;     % 尺度因子，即散射场振幅常系数与平面波常系数之间的比值
M_size = 401;   % 用于计算的矩阵的大小
for ii = 1:length(phi)      % 此循环用以生成不同phi值对应的散射条纹矩阵，并根据矩阵计算响应的psi值
    [Ei,Es,F,I] = wave_generate(lambda,n,kapa,theta,phi(ii),scale_factor,M_size);
    temp = I(ceil(M_size/2),ceil(M_size/2):end);    % Just take half of the middle section
    f = find_fft_bf(temp,1);
    f_test = @(para,xdata) cos(2*pi*para(1)*xdata + para(2));
    x = 0:(ceil(M_size/2)-1);
    y = temp./abs(hilbert(temp));
    lb = [0,0];ub = [2*pi,2*pi];    % set low and up bound for para
    opt = optimoptions('lsqcurvefit','Display','off');
    for jj = 1:20
        [p(jj).parameter,resnorm(jj)] = lsqcurvefit(f_test,[f,jj*pi/10],x(1:end-50),y(1:end-50),lb,ub,opt);
    end
    [~,loc] = min(resnorm);     
    psi(ii) = p(loc).parameter(2);
end
plot(phi,psi)

%% 解调法在两换不同交叉情况下的恢复结果
theta = 0:1:90;
n = [1.51,1.33];
lambda = 680;
kapa = 5;
phi = 0.1*pi;
scale_factor = 0.1;
M_size = 401;
mask = zeros(401);
mask(171:231,171:231) = 1;

savepath = uigetdir();
savepath = fullfile(savepath,'20211025_环离合与入射角的关系');
mkdir(savepath);
gifroute = fullfile(savepath,'不同较差情况下解调法恢复效果.gif');
for ii = 1:length(theta)
    [Ei,Es,F,I] = wave_generate(lambda,n,kapa,theta(ii),phi,scale_factor,M_size);
    temp = fftshift(fft2(exp(i*angle(Es.*conj(Ei))).*I));

    imagesc(abs(ifft2(ifftshift(mask.*temp))))
    axis off
    axis square
%     set(0,'defaultfigurecolor','w') 
    title(['incidence theta:',num2str(theta(ii)),'^o'],'fontsize',15,'fontweight','bold');
    F = getframe(gcf);
    I = frame2im(F);
    [I,map] = rgb2ind(I,256);
    if ii == 1
        imwrite(I,map,gifroute,'gif','Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map,gifroute,'gif','DelayTime',0.1,'WriteMode','append');
    end
    pause(0.01)
end

%%

tiffpath = 'G:\work\ScaterFeild\实验数据\20220227AuNPs_Ag_NWs_angle\20220227_2\TIFF\A4_NPs';
tiffs = dir(fullfile(tiffpath,'*.tiff'));
NoiTiffs = zeros(3,480,640);
for ii = 1:3
    NoiTiffs(ii,:,:) = double(imread(fullfile(tiffpath,tiffs(ii).name)));
end
[BG,ref] = eliminateBg(NoiTiffs,1e-3,200);
I = squeeze(BG(2,1:479,1:479));
figure
imagesc(I)
F = fftshift(fft2(I));
imagesc(log(abs(F)));
[mask,peaks,deg] = GetFourierMask(log(abs(F)),5,0,-1);


%% 载波
center = (size(I,1)+1)/2;   %  计算频域中心
ki = 2*pi*sqrt((peaks(1,1)-center)^2+(peaks(1,2)-center)^2)*1/size(I,1);    
% 通过圆心位置计算平面波的波矢，注意，我们以实空间的空间分辨率为1，那么k空间的频率频率分辨率就是1/size(I,1)了
ks = 2*pi*peaks(1,3)*1/size(I,1);
theta = -atan((peaks(2,2)-peaks(1,2))/(peaks(2,1)-peaks(1,1))); % 通过两个圆心的位置计算平面波传播方向
% theta = (360-deg)*pi/180;
Ec = zeros(size(I));
for ii = 1:size(I,1)
    for jj = 1:size(I,2)
        r = sqrt((ii-center)^2+(jj-center)^2);
        Ec(ii,jj) = exp(-i*(-ki*(ii-center)*sin(theta) + ki*(jj-center)*cos(theta)));
        Ec(ii,jj) = Ec(ii,jj)*exp(i*ks*r);
    end
end
figure
imagesc(angle(Ec))
figure
imagesc(log(abs(fftshift(fft2(Ec)))))

I_d = I.*Ec;
F_d = fftshift(fft2(I_d));
figure
imagesc((abs(F_d)))
    
mask = zeros(size(I));
mask(220:260,220:260) = 1;
figure
imagesc(abs(ifft2(ifftshift(mask.*F_d))))

