function bound=Training_LackOfRubber(pathstr)
% pathstr: the path name of the training set

JScoeff=[];
LengthDiff=[];
MeanGrad=[];
Sigm=[];
    for i=1:length(pathstr)
        % step 1: read and display the image 
        Img = double(imread(pathstr{i}));
        Img   = impreprocess(Img, 0.5, 10);
     
%         JS = check_excessofrubberThld(Img);
%         JScoeff = [JScoeff;JS];
        
        LD = check_excessofrubberThld2(Img);  
        LengthDiff = [LengthDiff;LD];
        
        MG = check_lackofrubberThld3(Img); %local mean gradient
        MeanGrad = [MeanGrad;MG];
        
        Sigm_value = check_lackofrubberLocalVari(Img); %local mean gradient
        Sigm = [Sigm;Sigm_value]; 
    end
    %     Up_JS=min(JScoeff); %low bound
    %     Low_LD=max(LengthDiff);%up bound
    %     Low_MG=max(MeanGrad);%up bound
%  	bound=[min(JScoeff) max(LengthDiff) max(MeanGrad) max(Sigm)];
bound = [max(LengthDiff) max(MeanGrad) max(Sigm)];
end