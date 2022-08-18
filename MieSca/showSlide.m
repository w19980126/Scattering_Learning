function [f]=showSlide(varargin)% showSlide(data,slidename,picnum,elsecommand,dimensionality)
    if isempty(varargin)
        msgbox({'Invalid Value',...
            'showSlide(data,slidename,picnum,elsecommand)'},...
            'Error','error');
        return;
    end
    
    data=varargin{1};
    varargin=varargin(2:end);
    slidename=('');
    picnum=1;
    elsecommand='';
    dimensionality=3;
    
    while length(varargin)>=2
        prop=varargin{1};
        val=varargin{2};
        varargin=varargin(3:end);
        eval([prop '=val;']);     
    end
    
    f=figure;
    drawModelEEMs(picnum,data,slidename,elsecommand,dimensionality);
    num=size(data,3);
    uicontrol('units','normalized','Style','slider','pos',[0 0 1 .05],...
        'min',1,'max',num,'value',1,...
        'sliderstep',[1/num,1/num],...
        'callback',@(obj,x,event)drawModelEEMs(round(get(obj,'value')),data,slidename,elsecommand,dimensionality));
end

function drawModelEEMs(ii,data,slidename,elsecommand,dimensionality)
    if dimensionality==1
        temp(:,:)=data(ii,:,:);
    end
    if dimensionality==2
        temp(:,:)=data(:,ii,:);
    end
    if dimensionality==3
        temp(:,:)=data(:,:,ii);
    end
    imagesc(temp);
    colormap(gray)
    title(['Slide ',num2str(ii),'-',slidename]);
    eval(elsecommand);
end

