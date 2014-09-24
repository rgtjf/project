function length_diff = check_excessofrubberThld2(Img)
% % segment the torus using threshold method 
display_fig = 0;
K = fspecial('gaussian',3,1);
Img = imfilter(Img,K,'same','symmetric');
[m,n,p] = size(Img);

g = double(rgb2gray(uint8(Img)));
IND = (g<mean(g(:))-50);
[uu,mask] = RemoveSmallComponent(1-double(IND), 0.5, 50);
IND = 1-uu;
% BWedge = edge(IND,'canny');% ʹ��canny edge�Ļ�,�е�edge����IND=1�ϣ��е���IND=0��
% ����ķ�������������edge�㶼��IND=1��

se = strel('disk',1);        
erodedBW = imerode(IND,se);
BWedge = IND-erodedBW;
% figure,imshow(BWedge)


BW = BWedge;BW2 = BW;
[row,col] = find(BW==1);
for kk=1:length(row)
aa = IND(row(kk),col(kk));
end
% figure,imshow(BW)
% % estimate the center of torus
[row,col] = find(BW==1);
y_center = mean(row);
x_center = mean(col);
if display_fig ==1;
figure,imshow(IND)
hold on; plot(x_center,y_center,'ro', 'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',5);hold off
end

% % find the inner and outer boundary add neibour by neibour
% ԭ���ķ����õ��ı߽������Ǵ�һ����ʼ�㿪ʼ����������ͬʱ����������x���겢������
% ��BW�г�����������
outerboundary_r = []; outerboundary_c = [];
[r,c] = find(BW,1,'first'); 
outerboundary_r = [outerboundary_r,r];
outerboundary_c = [outerboundary_c,c];
iter = 1;
outerboundaryBW = zeros(m,n);
outerboundaryBW(r,c)=1;
while iter<=length(outerboundary_r)
    iter = iter+1;  
    BW(r,c)=0; 
    outerboundaryBW(r,c)=1;    
    win = BW(max(r-1,1):min(r+1,m),max(c-1,1):min(c+1,n));
    [rr,cc] = find(win); 
    if ~isempty(rr)
       r = r+rr-2; c = c+cc-2;
       for l=1:length(r)
           BW(r(l),c(l))=0;  
           outerboundaryBW(r(l),c(l))=1;
       end
    outerboundary_r = [outerboundary_r,r'];
    outerboundary_c = [outerboundary_c,c'];   
    end
    if iter<=length(outerboundary_r)
    r = outerboundary_r(iter);
    c = outerboundary_c(iter);
    end
   
end

innerboundary_r1 = []; innerboundary_c1 = [];
innerboundary_r2 = []; innerboundary_c2 = [];

% [r1,c1] = find(BW,1,'first');
% [r2,c2] = find(BW(r1,:),1,'last');r2 = r1;
% [r3,c3] = find(BW,1,'last');
% [r4,c4] = find(BW(r3,:),1,'first');r4 = r3;
% figure; imshow(IND,[]);
% hold on; plot([c1,c2,c3,c4]',[r1,r2,r3,r4]','ro', 'MarkerEdgeColor','r',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',5);


[r,c] = find(BW,1,'first');  %�����кܶ��б�Ϊc�ĵ㣬�ҵ�������һ����Ϊ�����������
r = find(BW(:,c),1,'last');
BWU = zeros(size(BW));
BWU(1:r,:)=BW(1:r,:);% upper half
BWD = zeros(size(BW));
BWD(r+1:end,:)=BW(r+1:end,:);% bottom half
% figure,imshow([BWU BWD]);


[r,c] = find(BWU,1,'first');  %�����кܶ��б�Ϊc�ĵ㣬�ҵ�������һ����Ϊ�����������
r = find(BWU(:,c),1,'last');
BWU(r,c)=0;
innerboundary_r1 = [innerboundary_r1,r];
innerboundary_c1 = [innerboundary_c1,c];
iter = 1;
innerboundaryBW = zeros(m,n);
innerboundaryBW(r,c)=1;
while iter<=length(innerboundary_r1)
    iter = iter+1;
    BWU(r,c)=0;
    innerboundaryBW(r,c)=1;
    win = BWU(r-1:r+1,c-1:c+1);
    [rr,cc] = find(win);
    if ~isempty(rr)
       r = r+rr-2; c = c+cc-2;
       for l=1:length(r)
           BWU(r(l),c(l))=0;
           innerboundaryBW(r(l),c(l))=1;
       end
    innerboundary_r1 = [innerboundary_r1,r'];
    innerboundary_c1 = [innerboundary_c1,c'];
    end
    if iter<=length(innerboundary_r1)
    r = innerboundary_r1(iter);
    c = innerboundary_c1(iter);
    end

end
% ��������,��Ϊ��һ���ǶԳƵ�Բ�����԰�r1�ֳ���������Ļ�δ�������б굥�������Է�
% ����Ƭ������ߵĵ�����ұߵĵ�����б꽫���߷�Ϊ��С��
% innerboundary_rr =[innerboundary_r(innerboundary_r>=r1),innerboundary_r(innerboundary_r<r1)];
% innerboundary_cc =[innerboundary_c(innerboundary_r>=r1),innerboundary_c(innerboundary_r<r1)]

[r,c] = find(BWD,1,'first');  BWD(r,c)=0; 
innerboundary_r2 = [innerboundary_r2,r];
innerboundary_c2 = [innerboundary_c2,c];
iter = 1;
innerboundaryBW = zeros(m,n);
innerboundaryBW(r,c)=1;
while iter<=length(innerboundary_r2)
    iter = iter+1;  
    BWD(r,c)=0; 
    innerboundaryBW(r,c)=1;
    win = BWD(r-1:r+1,c-1:c+1);
    [rr,cc] = find(win); 
    if ~isempty(rr)
       r = r+rr-2; c = c+cc-2;
       for l=1:length(r)
           BWD(r(l),c(l))=0;  
           innerboundaryBW(r(l),c(l))=1;
       end
    innerboundary_r2 = [innerboundary_r2,r'];
    innerboundary_c2 = [innerboundary_c2,c'];   
    end
    if iter<=length(innerboundary_r2)
    r = innerboundary_r2(iter);
    c = innerboundary_c2(iter);
    end
   
end
innerboundary_r = [innerboundary_r1, fliplr(innerboundary_r2)];
innerboundary_c = [innerboundary_c1, fliplr(innerboundary_c2)];
if display_fig ==1;
figure; imshow(IND,[]);
for k=1:100:length(innerboundary_r)
hold on; plot(innerboundary_c(1:k),innerboundary_r(1:k),'ro', 'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',5);
            drawnow
end
end

            
            
            
% ������Ȧ���ⷨ��,��k������ڵ���Ϣ��PCAѰ���з���
hws = 30;
innerboundary_rr = padarray(innerboundary_r,[0 hws],'circular');
innerboundary_cc = padarray(innerboundary_c,[0,hws],'circular');

for i=1:length(innerboundary_r)
   rr = innerboundary_rr(i:i+2*hws)';
   cc = innerboundary_cc(i:i+2*hws)';
   X = cat(2,cc,rr);
   coeff = princomp(X,2);
  % coeff=coeff.A;
   nx(i)=coeff(1,2);
   ny(i)=-coeff(2,2);
end
k=length(innerboundary_c);
% hold on; quiver(innerboundary_c(1:k),innerboundary_r(1:k),nx(1:k),ny(1:k))
% axis equal
% % ע������������ķ����Ƕ�Ӧx���ң�y���ϵı�׼����ϵ�ģ�ֻ����ʾ��ͼ����ʱ���Զ���ת��
% % ����ͼ����Ӧ��Ϊ(nx,-ny)����;


% % ���㴩Խ�����ֱ������Ȧ��ֵmask BW2�Ľ��ߵĳ���
% ���ߵĲ�������
if isempty(innerboundary_r)% û���ڱ߽磬˵������Ȧ��ȱ��
    length_diff = 1000; 
else 
% line_length=func_line_intersection(IND,innerboundary_r,innerboundary_c,nx,ny);
line_length=func_line_intersection_R2(IND,innerboundary_r,innerboundary_c,nx,ny);

length_diff = max(line_length)-min(line_length);
end

% figure,plot(line_length)

end






% % ����ķ�������ķ���׼
% for i=1:length(innerboundary_r)
%    Tx = Dx(innerboundary_r); Ty = Dx(innerboundary_c);
%    T = sqrt(Tx.^2+Ty.^2);
%    T(T<1e-4) = 1e-4;
%    Tx = Tx./T; Ty = Ty./T;
% end
% hold on; quiver(innerboundary_c,innerboundary_r,Tx,Ty)

% 
% function d = Dx(u)
% [rows] = length(u); 
% d = zeros(1,rows);
% d(1,2:rows-1) = (u(1,3:rows)-u(1,1:rows-2))/2;
% return  

% [m,n] = size(IND);
% i = 1;
% x0 = innerboundary_c(i);
% y0 = innerboundary_r(i);
% for x=1:n
%     y = y0+tan(-ny(i)./nx(i))*(x-x0);
% end
% 
% 
% 
% t = 1;
% for t=1:100 
% x = x0+t*nx(i);
% y = y0+t*(-ny(i));
% if BW2(fix(x),fix(y))~=1
%    x = x0-t*nx(i);
%    y = y0-t*(-ny(i)); 
% end
% if BW2(fix(x),fix(y))==1
%    t = t+1;
% else
%     break;
% end
% end
% 
% 
% 
