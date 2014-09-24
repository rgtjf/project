% % 第二版的改进之处：原来统计像素个数，对于水平方向或者竖直方向的线来说是对的，
% % 但是对于斜线是不公平的，改为按照欧氏距离算
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

% 用y=kx+b斜率较大的线偏差很大
% 用参数方程
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
intersect = (line+IND==2); % 直线与torus交界线
intersectBW = intersect;
% figure,imshow(intersect);axis equal
% intersect包含x0，y0的联通分支即为所求
% 以x0,y0为起点，使用邻域搜索一个联通分支
r = y0; c = x0;
line_r = r;
line_c = c;
iter = 1;
lineBW = zeros(m,n);%记录联通分支，属于为1
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
% 求交线的长度，任意交线上两点算距离取最大值
[rr,cc] = find(lineBW==1);
for ii = 1:sum(lineBW(:))
    for jj = 1:sum(lineBW(:))
    temp_length(ii,jj) = sqrt((rr(ii)-rr(jj)).^2+(cc(ii)-cc(jj)).^2);
    end
end
intersect_line_length(i) = max(temp_length(:));
clear rr cc line line_r line_c x y intersectBW temp_length
end
