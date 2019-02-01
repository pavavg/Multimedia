function frameType = SSC(frameT, nextFrameT, prevFrameType)

b = [0.7548 -0.7548];
a = [1 -0.5095];
x1_filtered = filter(b,a,nextFrameT(:,1));
x2_filtered = filter(b,a,nextFrameT(:,2));

nextFrameType1 = 'OLS' ;
nextFrameType2 = 'OLS' ;

index = 449;
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

if prevFrameType == 'OLS'
    
    if nextFrameType1 == 'ESH'
        frameType1 = 'LSS' ;
    else 
        frameType1 = 'OLS' ;
    end
    
    if nextFrameType2 == 'ESH'
        frameType2 = 'LSS' ;
    else 
        frameType2 = 'OLS' ;
    end
    
elseif prevFrameType == 'ESH'
    
    if nextFrameType1 == 'ESH'
        frameType1 = 'ESH' ;
    else 
        frameType1 = 'LSS' ;
    end
    
    if nextFrameType2 == 'ESH'
        frameType2 = 'ESH' ;
    else 
        frameType2 = 'LSS' ;
    end
    
elseif prevFrameType == 'LSS'
    
    frameType1 = 'ESH' ;
    frameType2 = 'ESH' ;
    
else
    
    frameType1 = 'OLS' ;
    frameType2 = 'OLS' ;
    
end



if frameType1 == 'OLS'
    
    if frameType2 == 'OLS'
        frameType = 'OLS' ;
    elseif frameType2 == 'LSS'
        frameType = 'LSS' ;
    elseif frameType2 == 'ESH'
        frameType = 'ESH' ;
    else
        frameType = 'LSS' ; 
    end
    
elseif frameType1 == 'LSS'
    
    if frameType2 == 'OLS'
        frameType = 'LSS' ;
    elseif frameType2 == 'LSS'
        frameType = 'LSS' ;
    elseif frameType2 == 'ESH'
        frameType = 'ESH' ;
    else
        frameType = 'ESH' ; 
    end
    
elseif frameType1 == 'ESH'
    
    frameType = 'ESH' ;
    
else
    
    if frameType2 == 'OLS'
        frameType = 'LSS' ;
    elseif frameType2 == 'LSS'
        frameType = 'ESH' ;
    elseif frameType2 == 'ESH'
        frameType = 'ESH' ;
    else
        frameType = 'LSS' ; 
    end
    

end

