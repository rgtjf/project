function MG = check_lackofrubberThld3(Img)
% 将橡胶圈所在区域按圆心旋转一定角度分成若干份，计算每一份的
% 平均梯度，然后取最大值
display_fig = 0;
% % % % image preprocessing
K = fspecial('gaussian',3,1);
Img = imfilter(Img,K,'same','symmetric');


% % % % Segment the torus using threshold method
[m,n,p] = size(Img);
g = double(rgb2gray(uint8(Img)));
BW = (g<mean(g(:))-50);
[uu,mask] = RemoveSmallComponent(1-double(BW), 0.5, 50);
IND = 1-uu;
% BWedge = edge(IND,'canny');% 使用canny edge的话,有的edge点在IND=1上，有的在IND=0上
% 下面的方法可以让所有edge点都在IND=1上
se = strel('disk',6);        
BW = imerode(IND,se);
if display_fig==1
figure,imshow(uint8(Img));
hold on; contour(BW,[0.5,0.5],'r')
end
% 找torus中心，将[0，2PI]等分成若干份
[y_center, x_center] = find_center(BW);
[y,x] = meshgrid(1:m,1:n);
y = y-y_center;
x = x-x_center;
theta = atan2(x,y)+pi;
% figure,imshow(theta',[])
n = 10; % subregions_number
theta0 = 2*pi/n;
if display_fig==1
figure
for i = 1:n
    subregion_BW(:,:,i) = BW.*(theta'>(i-1)*theta0& theta'<i*theta0);
    imshow(subregion_BW(:,:,i),[]);drawnow
% calculate the gradient on each subregion
grad(:,:,i) = sqrt(sum(Dx(Img).^2+Dy(Img).^2,3)).*subregion_BW(:,:,i);
meangrad(i) = mean2(grad(:,:,i));
end
else
  for i = 1:n
    subregion_BW(:,:,i) = BW.*(theta'>(i-1)*theta0& theta'<i*theta0);
    grad(:,:,i) = sqrt(sum(Dx(Img).^2+Dy(Img).^2,3)).*subregion_BW(:,:,i);    
    meangrad(i) = mean2(grad(:,:,i));
  end
end
MG = max(meangrad);

function d = Dx(u)
[rows,cols,p] = size(u); 
d = zeros(rows,cols,p);
d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
d(:,1,:) = u(:,1,:)-u(:,cols,:);
return

function d = Dy(u)
[rows,cols,p] = size(u); 
d = zeros(rows,cols,p);
d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
d(1,:,:) = u(1,:,:)-u(rows,:,:);
return


