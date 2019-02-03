function AACSeq1 = AACoder1(fNameIn)
y = audioread(fNameIn);
frameNumber = floor(size(y,1) /1024 -2);

AACSeq1 = struct([]) ;
index = 1;
prevFrameType = 'OLS' ;
for i=1:frameNumber
    frameT = y(index:index+2047 , :) ;
    nextFrameT = y(index+1024: index+3071, :) ;
    
    AACSeq1(i).frameType = SSC(frameT, nextFrameT, prevFrameType) ;
    prevFrameType = AACSeq1(i).frameType;
    
    AACSeq1(i).winType = 'KBD' ;
    
    frameF = filterbank(frameT, AACSeq1(i).frameType, AACSeq1(i).winType);
    
    AACSeq1(i).chl.frameF = frameF(:,1) ;
    AACSeq1(i).chr.frameF = frameF(:,2) ;
    
    index = index +1024;
end

