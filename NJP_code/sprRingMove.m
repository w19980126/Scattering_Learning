%% Step 1: 如果两圆心不在垂直线上，则需要旋转插值，使其呈8字

%% Step 2: 估计两个图像上环移动的距离，即计算dist

%% Step 3: 垂直方向移动图像拆分
%示例数据生成
siz=401;
im=peaks(siz);
IM=im+rot90(im,2);
dist=10*1;
temp=eye(siz);
moveMat=[temp(dist+1:end,:);temp(1:dist,:)];
IM2=moveMat*im+rot90(moveMat*im,2);
IM2=im*moveMat+rot90(im*moveMat,2);

%拆分
res = separateIM(IM,IM2,dist,'r');

%画图
subplot(3,2,1)
imagesc(IM);title('旋转合并图像1');
subplot(3,2,2)
imagesc(IM2);title('旋转合并图像2');
subplot(3,2,3)
imagesc(im);title('原始图像');
subplot(3,2,4)
imagesc(res);title('还原图像');
subplot(3,2,5)
imagesc((im-res)./im);title('偏差');

IM = IM1';
IM2 = IM2';
dist = 9;
%画图
subplot(3,2,1)
imagesc(abs(IM));title('旋转合并图像1');
subplot(3,2,2)
imagesc(abs(IM2));title('旋转合并图像2');
subplot(3,2,3)
imagesc(abs(IM));title('原始图像');
subplot(3,2,4)
imagesc(abs(res));title('还原图像');
subplot(3,2,5)
imagesc(abs(im-res)./abs(im));title('偏差');










