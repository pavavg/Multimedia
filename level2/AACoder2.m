function AACSeq2 = AACoder2(fNameIn)

y = audioread(fNameIn); %Read Audio
frameNumber = floor(size(y,1) /1024 -2); %Frame number to be splitted

AACSeq2 = struct([]) ;
index = 1;
prevFrameType = 'OLS' ;

for i=1:frameNumber
    frameT = y(index:index+2047 , :) ;
    nextFrameT = y(index+1024: index+3071, :) ;
    
    AACSeq2(i).frameType = SSC(frameT, nextFrameT, prevFrameType) ;
    prevFrameType = AACSeq2(i).frameType;
    
    AACSeq2(i).winType = 'SIN' ;
    
    frameF = filterbank(frameT, AACSeq2(i).frameType, AACSeq2(i).winType);
    
    [AACSeq2(i).chl.frameF ,AACSeq2(i).chl.TNScoeffs] = TNS(frameF(:,1,:), AACSeq2(i).frameType) ; 
    [AACSeq1(i).chr.frameF ,AACSeq2(i).chr.TNScoeffs] = TNS(frameF(:,2,:), AACSeq2(i).frameType) ;
    
    index = index +1024;
end

%Last frame
frameT = y(index:index+2047 , :) ;
AACSeq2(frameNumber+1).frameType = AACSeq2(frameNumber).frameType ;
AACSeq2(frameNumber+1).winType = 'SIN' ;
frameF = filterbank(frameT, AACSeq2(frameNumber+1).frameType, AACSeq2(frameNumber+1).winType);
[AACSeq2(frameNumber+1).chl.frameF ,AACSeq2(frameNumber+1).chl.TNScoeffs] = TNS(frameF(:,1,:), AACSeq2(frameNumber+1).frameType) ; 
[AACSeq1(frameNumber+1).chr.frameF ,AACSeq2(frameNumber+1).chr.TNScoeffs] = TNS(frameF(:,2,:), AACSeq2(frameNumber+1).frameType) ;
    
end
