function frameType = SSC(frameT, nextFrameT, prevFrameType)

b = [0.7548 -0.7548];
a = [1 -0.5095];
x1_filtered = filter(b,a,nextFrameT(:,1));
x2_filtered = filter(b,a,nextFrameT(:,2));

nextFrameType1 = 'OLS' ;
nextFrameType2 = 'OLS' ;

index = 449+128;
s_sq1 = zeros(8,1);
s_sq2 = zeros(8,1);

ds_sq1 = zeros(8,1);
ds_sq2 = zeros(8,1);
for i=1:8
    subframe1 = x1_filtered(index:index+127);
    subframe2 = x2_filtered(index:index+127);
    index = index + 128;
    s_sq1(i) = sum( subframe1.^2);
    s_sq2(i) = sum( subframe2.^2);
    
    if i>1 && s_sq1(i) > 0.001
        ds_sq1(i) = s_sq1(i) / mean( s_sq1(1:i-1) );
        
        if ds_sq1(i) > 10
            nextFrameType1 = 'ESH' ;
            break;
        end 
    end
    
    if i>1 && s_sq2(i) > 0.001
        ds_sq2(i) = s_sq2(i) / mean( s_sq2(1:i-1) );
        
        if ds_sq2(i) > 10
            nextFrameType2 = 'ESH' ;
            break;
        end
    end
    
end

if strcmp(prevFrameType,'OLS')
    
    if strcmp(nextFrameType1 ,'ESH')
        frameType1 = 'LSS' ;
    else 
        frameType1 = 'OLS' ;
    end
    
    if strcmp(nextFrameType2 ,'ESH')
        frameType2 = 'LSS' ;
    else 
        frameType2 = 'OLS' ;
    end
    
elseif strcmp(prevFrameType, 'ESH')
    
    if strcmp(nextFrameType1 ,'ESH')
        frameType1 = 'ESH' ;
    else 
        frameType1 = 'LPS' ;
    end
    
    if strcmp(nextFrameType2 ,'ESH')
        frameType2 = 'ESH' ;
    else 
        frameType2 = 'LPS' ;
    end
    
elseif strcmp(prevFrameType ,'LSS')
    
    frameType1 = 'ESH' ;
    frameType2 = 'ESH' ;
    
else
    
    frameType1 = 'OLS' ;
    frameType2 = 'OLS' ;
    
end



if strcmp(frameType1 ,'OLS')
    
    if strcmp(frameType2 ,'OLS')
        frameType = 'OLS' ;
    elseif strcmp(frameType2 ,'LSS')
        frameType = 'LSS' ;
    elseif strcmp(frameType2 ,'ESH')
        frameType = 'ESH' ;
    else
        frameType = 'LPS' ; 
    end
    
elseif strcmp(frameType1 ,'LSS')
    
    if strcmp(frameType2 ,'OLS')
        frameType = 'LSS' ;
    elseif strcmp(frameType2 ,'LSS')
        frameType = 'LSS' ;
    elseif strcmp(frameType2 ,'ESH')
        frameType = 'ESH' ;
    else
        frameType = 'ESH' ; 
    end
    
elseif strcmp(frameType1 ,'ESH')
    
    frameType = 'ESH' ;
    
else
    
    if strcmp(frameType2 ,'OLS')
        frameType = 'LPS' ;
    elseif strcmp(frameType2 ,'LSS')
        frameType = 'ESH' ;
    elseif strcmp(frameType2 ,'ESH')
        frameType = 'ESH' ;
    else
        frameType = 'LPS' ; 
    end
    

end

