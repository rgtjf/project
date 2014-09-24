% % �ڶ���ĸĽ�֮����ԭ��ͳ�����ظ���������ˮƽ���������ֱ���������˵�ǶԵģ�
% % ���Ƕ���б���ǲ���ƽ�ģ���Ϊ����ŷ�Ͼ�����
function intersect_line_length = func_line_intersection_R2(IND,innerboundary_r,innerboundary_c,nx,ny)
for i=1:length(innerboundary_c)
IND(innerboundary_r(i),innerboundary_c(i))=1;    
end
% figure,imshow(IND)
% k=1;
% hold on; plot(innerboundary_c(k),innerboundary_r(k),'ro', 'MarkerEdgeColor','r',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',5);
% hold on; quiver(innerboundary_c(k),innerboundary_r(k),nx(k),ny(k),'LineWidth',5)
% axis equal

BW2 = IND;
[m,n] = size(IND);

% for i=316
for i=1:length(innerboundary_c)
x0 = innerboundary_c(i);
y0 = innerboundary_r(i);

% ��y=kx+bб�ʽϴ����ƫ��ܴ�
% �ò�������
t = -100:1:100;
for k=1:length(t)
x(k) = x0+t(k)*nx(i);
y(k) = y0+t(k)*ny(i);
end
xx = x(x>=1&x<n&y>=1&y<m);
yy = y(x>=1&x<n&y>=1&y<m);
x = xx;
y = yy;

% hold on; plot(x,y,'g','LineWidth',5)
line = zeros(m,n);
for kk = 1:length(y)    
line(floor(y(kk)),floor(x(kk))) = 1;    
end
intersect = (line+IND==2); % ֱ����torus������
intersectBW = intersect;
% figure,imshow(intersect);axis equal
% intersect����x0��y0����ͨ��֧��Ϊ����
% ��x0,y0Ϊ��㣬ʹ����������һ����ͨ��֧
r = y0; c = x0;
line_r = r;
line_c = c;
iter = 1;
lineBW = zeros(m,n);%��¼��ͨ��֧������Ϊ1
lineBW(r,c)=1;
while iter<=length(line_r)
    iter = iter+1;  
    intersectBW(r,c)=0; 
    lineBW(r,c)=1;
    win = intersectBW(max(r-1,1):min(r+1,m),max(c-1,1):min(c+1,n));
    [rr,cc] = find(win); 
    if ~isempty(rr)
       r = r+rr-2; c = c+cc-2;
       for l=1:length(r)
           intersectBW(r(l),c(l))=0;  
           lineBW(r(l),c(l))=1;
       end
    line_r = [line_r,r'];
    line_c = [line_c,c'];   
    end
    if iter<=length(line_r)
    r = line_r(iter);
    c = line_c(iter);
    end
end
% ���ߵĳ��ȣ����⽻�������������ȡ���ֵ
[rr,cc] = find(lineBW==1);
for ii = 1:sum(lineBW(:))
    for jj = 1:sum(lineBW(:))
    temp_length(ii,jj) = sqrt((rr(ii)-rr(jj)).^2+(cc(ii)-cc(jj)).^2);
    end
end
intersect_line_length(i) = max(temp_length(:));
clear rr cc line line_r line_c x y intersectBW temp_length
end
