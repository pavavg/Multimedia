function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
load('TableB219.mat')
forceCodebook = 12;
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
    SMRl = psycho(frameT(:,1), AACSeq3(i).frameType ,frameTprev1(:,1), frameTprev2(:,1) );
    SMRr = psycho(frameT(:,2), AACSeq3(i).frameType ,frameTprev1(:,2), frameTprev2(:,2) );
    
    if ~strcmp(AACSeq3(i).frameType, 'ESH')
        noBands = size(B219a, 1) ;
        Pl = zeros(noBands,1);
        Pr = zeros(noBands,1);
        for b = 1:noBands
            start = B219a(b, 2)+1;
            finish = B219a(b, 3) +1;
            Pl(b)  = sum(AACSeq3(i).chl.frameF(start:finish,: ).^2) ;
            Pr(b)  = sum(AACSeq3(i).chr.frameF(start:finish,: ).^2) ;
        end
        AACSeq3(i).chl.T = Pl ./SMRl ;
        AACSeq3(i).chr.T = Pr ./SMRr ;
    else
        noBands = size(B219b, 1) ;
        Pl = zeros(noBands,8);
        Pr = zeros(noBands,8);
        
        for f =1:8
            for b = 1:noBands
                start = B219b(b, 2)+1;
                finish = B219b(b, 3) +1;
                Pl(b, f)  = sum(AACSeq3(i).chl.frameF(start: finish ,f).^2) ;
                Pr(b, f)  = sum(AACSeq3(i).chr.frameF(start: finish ,f).^2) ;
            end
        end
        AACSeq3(i).chl.T = Pl ./SMRl ;
        AACSeq3(i).chr.T = Pr ./SMRr ;
        
    end
    
    [Sl, sfcL, AACSeq3(i).chl.G] = AACquantizer(frameF(:,1,:), AACSeq3(i).frameType, SMRl);
    [Sr, sfcR, AACSeq3(i).chr.G] = AACquantizer(frameF(:,2,:), AACSeq3(i).frameType, SMRr);
    if ~strcmp(AACSeq3(i).frameType, 'ESH')
        [AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook] = encodeHuff(Sl, loadLUT() );
        [AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook] = encodeHuff(Sr, loadLUT() );
        %[AACSeq3(i).chl.sfc, ~] = encodeHuff(sfcL, loadLUT(), forceCodebook);
        %[AACSeq3(i).chr.sfc,~] = encodeHuff(sfcR, loadLUT() , forceCodebook);
    else
        streamL = [];
        streamR = [];
        for f =1:8
            [tempStream, tempCodebookL] = encodeHuff(Sl(:,f), loadLUT() );
            streamL = strcat(streamL,tempStream) ;
            [tempStream, tempCodebookR] = encodeHuff(Sr(:,f), loadLUT() );
            streamR = strcat(streamR ,tempStream) ;
        end
        AACSeq3(i).chl.stream = streamL;
        AACSeq3(i).chr.stream = streamR;
        AACSeq3(i).chl.codebook = tempCodebookL;
        AACSeq3(i).chr.codebook = tempCodebookR;
    end
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

