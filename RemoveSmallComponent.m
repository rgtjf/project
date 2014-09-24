function [uu,mask] = RemoveSmallComponent(u, thr, mmsize)

graph=u/(max(u(:))-min(u(:)))>thr;
BW = logical(graph);
[L, num] = bwlabel(BW, 8);

ccsize=zeros(num,1);
for k=1:num
   [r, c] = find(bwlabel(BW)==k);
   %ccsize(k)=max(size(r));
   ccsize(k)=max(max(r)-min(r), max(c)-min(c));
end

mask=zeros(size(L));
for i=1:size(L,1)
  for j=1:size(L,2)
    mask(i,j)= double( (L(i,j)>0) &&( ccsize(L(i,j))>mmsize));
    
  end
end

uu=u.*mask;
end