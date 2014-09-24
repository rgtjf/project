function Img = impreprocess(Image,re_size,paddsize)
% image processing :
% 1. reduce the size
% 2. cut blank edges

% re_size  调整图像尺寸
% paddsize 扩展边界
Image=imresize(Image,re_size);
Max_Img=max(Image(:));
[w,h,c]=size(Image);
%  GImage=rgb2gray(Image./255);
Mean_Img=mean(Image(:));

while 1
    BreakLabel=0;
    if sum(sum(Image(:,1,:)>Mean_Img))>=c*w
        Image=Image(:,2:end,:);
        h=h-1;
    else
        BreakLabel=BreakLabel+1;
    end
    
    if sum(sum(Image(:,h,:)>Mean_Img))>=c*w
        Image=Image(:,1:end-1,:);
        h=h-1;
    else
        BreakLabel=BreakLabel+1;
    end
    
     if sum(sum(Image(1,:,:)>Mean_Img))>=c*h
        Image=Image(2:end,:,:);
        w=w-1;
    else
        BreakLabel=BreakLabel+1;
     end
     
     if sum(sum(Image(w,:,:)>Mean_Img))>=c*h
        Image=Image(1:end-1,:,:);
        w=w-1;
    else
        BreakLabel=BreakLabel+1;
     end
    
     if BreakLabel==4
         break
     end
%      imshow(Image./255)
end 
Img=ones(w+2*paddsize,h+2*paddsize,c)*double(Max_Img);
Img(paddsize+1:w+paddsize,paddsize+1:h+paddsize,:)=Image; 

end