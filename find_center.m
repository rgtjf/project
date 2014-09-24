function [y_center, x_center] = find_center(IND)
display_fig = 0;
BW = edge(IND,'canny');
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
