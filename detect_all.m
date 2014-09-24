function reslt=detect_all()

close all; clear; clc;

% add the source folder to the search path
path('src',path); 
path('src/lackrubber',path);
path('src/excessrubber',path);

% search all the files to check
pth = '20140627\合格';
[~,pathstr] = dirext(pth,1,'bmp');

bound=Training_LackOfRubber(pathstr);


pth = '20140627\毛边';
[~,pathstr] = dirext(pth,1,'bmp');
for i=1:length(pathstr)
    reslt(i)=detect_LackOfRubber(pathstr{i},bound); % 0 is right
end

pth = '20140627\缺料';
[~,pathstr] = dirext(pth,1,'bmp');
for j=1:length(pathstr)
    reslt(i+j)=detect_LackOfRubber(pathstr{j},bound); % 0 is right 1 is bad
end

 pth = '20140627\合格';
[~,pathstr] = dirext(pth,1,'bmp');
for k=1:length(pathstr)
reslt(i+j+k)=detect_LackOfRubber(pathstr{k},bound); % 0 is right 1 is bad
end

%precision_rate=sum([reslt(1:20) ~reslt(21:end)])/(i+j+k);
end