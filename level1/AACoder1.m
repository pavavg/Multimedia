function AACSeq1 = AACoder1(fNameIn)
y = audioread(fNameIn); %Read Audio
frameNumber = floor(size(y,1) /1024 -2); %Frame number to be splitted

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
    
    AACSeq1(i).chl.frameF = frameF(:,1,:) ;
    AACSeq1(i).chr.frameF = frameF(:,2,:) ;
    
    index = index +1024;
end

%Last frame
frameT = y(index:index+2047 , :) ;
AACSeq1(frameNumber+1).frameType = AACSeq1(frameNumber).frameType ;
AACSeq1(frameNumber+1).winType = 'KBD' ;
frameF = filterbank(frameT, AACSeq1(frameNumber+1).frameType, AACSeq1(frameNumber+1).winType);
AACSeq1(frameNumber+1).chl.frameF = frameF(:,1,:) ;
AACSeq1(frameNumber+1).chr.frameF = frameF(:,2,:) ;

end
