% 图像的傅里叶变换
% 2021.12.26

function FIm=PictureFourier(Im)
if length(size(Im))==2
    N=1;
elseif length(size(Im))==3
    [N,~,~]=size(Im);
end
if N==1
    fftI=fftshift(fft2(squeeze(Im)));
    siz=size(fftI);
    sizMax=max(siz);
    [X,Y]=meshgrid(1/sizMax:siz(2)/sizMax:siz(2),1/sizMax:siz(1)/sizMax:siz(1));
    FIm=interp2(log(abs(fftI)),X,Y);
else
    fftI1=fftshift(fft2(squeeze(Im(1,:,:))));
    siz=size(fftI1);
    sizMax=max(siz);
    [X,Y]=meshgrid(1/sizMax:siz(2)/sizMax:siz(2),1/sizMax:siz(1)/sizMax:siz(1));
    FIm1=interp2(log(abs(fftI1)),X,Y);
    [a,b]=size(FIm1);
    FIm=zeros(N,a,b);
    FIm(1,:,:)=FIm1;
    for i=2:N
        fftI=fftshift(fft2(squeeze(Im(i,:,:))));
        siz=size(fftI);
        sizMax=max(siz);
        [X,Y]=meshgrid(1/sizMax:siz(2)/sizMax:siz(2),1/sizMax:siz(1)/sizMax:siz(1));
        FIm(i,:,:)=interp2(log(abs(fftI)),X,Y);
    end
end
end