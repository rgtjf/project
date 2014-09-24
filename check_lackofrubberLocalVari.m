function MG = check_lackofrubberLocalVari(Img)
% ����Ȧ��������Բ����תһ���Ƕȷֳ�n�ݣ�����ÿһ�ݵ�
% ���Ȼ��ȡ���ֵ������n�ж���ȡֵ

% % % % image preprocessing
K = fspecial('gaussian',3,1);
Img = imfilter(Img,K,'same','symmetric');


% % % % Segment the torus using threshold method
[m,n,p] = size(Img);
g = double(rgb2gray(uint8(Img)));
BW = (g<mean(g(:))-50);
[uu,mask] = RemoveSmallComponent(1-double(BW), 0.5, 50);
IND = 1-uu;
% BWedge = edge(IND,'canny');% ʹ��canny edge�Ļ�,�е�edge����IND=1�ϣ��е���IND=0��
% ����ķ�������������edge�㶼��IND=1��
se = strel('disk',6);        
BW = imerode(IND,se);
% if display_fig==1
% figure,imshow(uint8(Img));
% hold on; contour(BW,[0.5,0.5],'r')
% end
% ��torus���ģ���[0��2PI]�ȷֳ����ɷ�
[y_center, x_center] = find_center(BW);
[y,x] = meshgrid(1:m,1:n);
y = y-y_center;
x = x-x_center;
theta = atan2(x,y)+pi;
% figure,imshow(theta',[])

kk=1;
for n=[16 20 24 30 36 54]
%n = 10; % subregions_number
theta0 = 2*pi/n;

 for i = 1:n
    subregion_BWi = BW.*(theta'>(i-1)*theta0& theta'<i*theta0);
    subregion_BWi3=cat(3,subregion_BWi,subregion_BWi,subregion_BWi);
   
    subregion_Img=Img(subregion_BWi3(:)>0);
    sigma(i)=var(subregion_Img);
    
%       subregion_Img = Img.*subregion_BWi3; 
%     mean_subregion = sum( subregion_Img(:) )./sum( subregion_BWi3(:) );
%     sigma_subregion = (subregion_Img(:)-mean_subregion).^2.*subregion_BWi3(:);
%     sigma(i)=sqrt(sum( sigma_subregion(:) ) ./ sum( subregion_BWi3(:) )); 
end
MG(kk) = max(sigma);kk=kk+1;
end


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


