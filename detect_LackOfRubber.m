function reslt=detect_LackOfRubber(Im_name,bound)
% pathstr: the path name of the training set
 
   % step 1: read and display the image 
    Img = double(imread(Im_name));
    Img   = impreprocess(Img, 0.5, 10);
   
%     JS = check_excessofrubberThld(Img);
%     if JS<bound(1)
%         reslt=1;  % if ==1 then this image is a bad one
%         return;
%     end
    
    LD = check_excessofrubberThld2(Img); %the length form 
    if LD>bound(1) 
        reslt=1;
        return;
    end
%     
%     MG = check_lackofrubberThld2(Img);
%     if MG>bound(2)
%          reslt=1;
%         return;
%     end
    
    Sigm=check_lackofrubberLocalVari(Img); 
    if sum( Sigm>bound(3:end) )>0
         reslt=1;
        return;
    end
    
    reslt=0;
end