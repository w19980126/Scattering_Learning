function [mask,peaks,deg]=GetFourierMask(preMask2,lineWidth,deg,method) % 以手动方式确定mask
% method:负数：多点确认半径；0：三点确认半径
% deg：预设的圆环偏离水平的角度，用角度制
% linewidth：圆环粗细
% preMasks：输入的图片
siz=size(preMask2);
sizMax=max(siz);
if method==0
    figure;imagesc(preMask2);
    title('Click 3 points on one circle');
    [x1,y1]=ginput(3);
    [x1,y1]=regulatepoint(x1,y1,preMask2,2);% 画点后在周围找到最合适的点
    title('Click 3 points on another circle');
    [x2,y2]=ginput(3);
    [x2,y2]=regulatepoint(x2,y2,preMask2,2);
    [c1,r1]=calc(x1,y1); % 三点确认圆心与半径，并画出来
    [c2,r2]=calc(x2,y2);
else
    preMask2(isnan(preMask2))=0;preMask2=to01(preMask2);
    figure;imagesc(preMask2);
    pause(10);
    while 1
        answer=questdlg('Continue?','Continue?','YES','NO','YES');
        switch answer
            case 'YES'
                break;
            case 'NO'
                pause(10);
        end
    end
    while 1
        prompt = {'level:','block length:'};
        dlgtitle = 'Input';
        definput = {'0.8','10'};
        dims = [1 20];
        answer=inputdlg(prompt,dlgtitle,dims,definput);
        answer=str2double(answer);
        level=answer(1);
        blocklength=answer(2);
        Block=ones(siz);
        Block(:,round(siz(1)/2-blocklength):round(siz(1)/2+blocklength))=0;
        BW=im2bw(preMask2,level);
        BW=BW.*Block;
        F=figure;imagesc(BW);
        pause(3);
        answer=questdlg('Continue?','Continue?','YES','NO','YES');
        switch answer
            case 'YES'
                close(F);
                break;
            case 'NO'
                close(F);
        end
    end
    [pointy,pointx]=find(BW==1);
    y1=pointy(pointx<round(siz(1)/2));
    x1=pointx(pointx<round(siz(1)/2));
    y2=pointy(pointx>round(siz(1)/2));
    x2=pointx(pointx>round(siz(1)/2));
    [c1,r1]=Multicalc(x1,y1);
    [c2,r2]=Multicalc(x2,y2);    
end
caltan=(c1(2)-c2(2))/(c1(1)-c2(1));
if c1(1)==c2(1) % 去除奇点
    adeg=90;
    if deg>180
        adeg=270;
    end
    deg=adeg;
elseif c1(2)==c2(2)
    adeg=0;
    if deg>90
        adeg=180;
    end
    deg=adeg;
else
    if caltan>0
        adeg=atand(caltan);
        if deg>135 && deg<225
            adeg=adeg+180;
        end
        deg=adeg;
    else
        adeg=360+atand(caltan);
        if deg>135 && deg<225
            adeg=adeg-180;
        end
        deg=adeg;
    end
end
r=(r1+r2)/2;
peaks=[c1,r;c2,r];
maskTemp=zeros(size(preMask2));
temp=insertShape(maskTemp,'circle',peaks,'LineWidth',lineWidth,'Color','white');
maskTemp=temp(:,:,1);
[X,Y]=meshgrid(1/siz(2):sizMax/siz(2):sizMax,1/siz(1):sizMax/siz(1):sizMax);
mask = interp2(maskTemp,X,Y);
mask(1,:)=0;
mask(:,1)=0;
end

function [cc,r]=Multicalc(x,y) % 多点确认圆心与半径，并画出来
n=length(x);
xx=x.*x;
yy=y.*y;
xy=x.*y;
A=[sum(x) sum(y) n ; sum(xy) sum(yy) sum(y) ; sum(xx) sum(xy) sum(x) ];
B=[-sum(xx+yy) ; -sum(xx.*y+yy.*y) ; -sum(xx.*x+xy.*y) ];
a=A\B;
xc=-a(1)/2;
yc=-a(2)/2;
cc=[xc yc];
r=sqrt((a(1)^2+a(2)^2)/4-a(3));
deg=0:360;           
rx=cc(1)+r*cosd(deg);
ry=cc(2)+r*sind(deg); 
hold on ;  
plot(x,y,'ko'); 
plot(rx,ry,'r');hold off;
pause(1);
end

function IM=to01(im)
immax=max(max(im));
immin=min(min(im));
IM=(im-immin)./(immax-immin);
end