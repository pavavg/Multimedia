function AACSeq3 = AACoder3(fNameIn, fnameAACoded)

y = audioread(fNameIn); %Read Audio
frameNumber = floor(size(y,1) /1024 -2); %Frame number to be splitted

AACSeq3 = struct([]) ;
index = 1;
prevFrameType = 'OLS' ;
frameTprev1 = zeros(2048,2);
frameTprev2 = zeros(2048,2);
for i=1:frameNumber
    frameT = y(index:index+2047 , :) ;
    nextFrameT = y(index+1024: index+3071, :) ;
    
    AACSeq3(i).frameType = SSC(frameT, nextFrameT, prevFrameType) ;
    prevFrameType = AACSeq3(i).frameType;
    
    AACSeq3(i).winType = 'SIN' ;
    
    frameF = filterbank(frameT, AACSeq3(i).frameType, AACSeq3(i).winType);
    if  strcmp(AACSeq3(i).frameType ,'ESH')
       % size(frameF)
    end
    [AACSeq3(i).chl.frameF ,AACSeq3(i).chl.TNScoeffs] = TNS(frameF(:,1,:), AACSeq3(i).frameType) ; 
    [AACSeq3(i).chr.frameF ,AACSeq3(i).chr.TNScoeffs] = TNS(frameF(:,2,:), AACSeq3(i).frameType) ;
    SMRr = psycho(frameT(:,1), AACSeq3(i).frameType ,frameTprev1(:,1), frameTprev2(:,1) );
    SMRl = psycho(frameT(:,2), AACSeq3(i).frameType ,frameTprev1(:,2), frameTprev2(:,2) );
    [S, sfc, G] = AACquantizer(frameF(:,1,:), AACSeq3(i).frameType, SMRr);
    
    frameTprev2 = frameTprev1;
    frameTprev1 = frameT;
    index = index +1024;
end

%Last frame
frameT = y(index:index+2047 , :) ;
AACSeq3(frameNumber+1).frameType = AACSeq3(frameNumber).frameType ;
AACSeq3(frameNumber+1).winType = 'SIN' ;
frameF = filterbank(frameT, AACSeq3(frameNumber+1).frameType, AACSeq3(frameNumber+1).winType);
[AACSeq3(frameNumber+1).chl.frameF ,AACSeq3(frameNumber+1).chl.TNScoeffs] = TNS(frameF(:,1,:), AACSeq3(frameNumber+1).frameType) ; 
[AACSeq3(frameNumber+1).chr.frameF ,AACSeq3(frameNumber+1).chr.TNScoeffs] = TNS(frameF(:,2,:), AACSeq3(frameNumber+1).frameType) ;
    
save(fnameAACoded, 'AACSeq3')
end

