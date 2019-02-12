function x = iAACoder3(AACSeq3, fNameOut)

frameNumber = size(AACSeq3 , 2);

channel_1 = [];
channel_2 = [];

sfcL = decodeHuff(AACSeq3(1).chl.sfc, 12, loadLUT());
sfcR = decodeHuff(AACSeq3(1).chr.sfc, 12, loadLUT());
S_L = decodeHuff(AACSeq3(1).chl.stream, AACSeq3(1).chl.codebook, loadLUT());
S_R = decodeHuff(AACSeq3(1).chr.stream, AACSeq3(1).chr.codebook, loadLUT());
frameF_L = iAACquantizer(S_L, [AACSeq3(1).chl.G ; sfcL'], AACSeq3(1).chl.G, AACSeq3(1).frameType);
frameF_R = iAACquantizer(S_R, [AACSeq3(1).chr.G ; sfcR'], AACSeq3(1).chr.G, AACSeq3(1).frameType);

frameF = [iTNS(frameF_L, AACSeq3(1).frameType, AACSeq3(1).chl.TNScoeffs) iTNS(frameF_R, AACSeq3(1).frameType, AACSeq3(1).chr.TNScoeffs)] ;

frameT = iFilterbank( frameF, AACSeq3(1).frameType, AACSeq3(1).winType );
prevFrameT = frameT;


for i=2:frameNumber
    
    if strcmp(AACSeq3(i).frameType ,'ESH')
        sfcL = decodeHuff(AACSeq3(i).chl.sfc, 12, loadLUT());
        sfcR = decodeHuff(AACSeq3(i).chr.sfc, 12, loadLUT());
        S_L = decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook, loadLUT());
        S_R = decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook, loadLUT());
        
        sfcL = reshape(sfcL,41,8);
        sfcR = reshape(sfcR,41,8);
        S_L = reshape(S_L,128,8);
        S_R = reshape(S_R,128,8);
        
        
        frameF_L = iAACquantizer(S_L, [AACSeq3(i).chl.G ; sfcL], AACSeq3(i).chl.G, AACSeq3(i).frameType);
        frameF_R = iAACquantizer(S_R, [AACSeq3(i).chr.G ; sfcR], AACSeq3(i).chr.G, AACSeq3(i).frameType);
        frameF = cat(3, iTNS(frameF_L, AACSeq3(i).frameType, AACSeq3(i).chl.TNScoeffs) , iTNS(frameF_R, AACSeq3(i).frameType, AACSeq3(i).chr.TNScoeffs)) ;
        
        frameF = permute(frameF, [1 3 2]);
    else
        sfcL = decodeHuff(AACSeq3(i).chl.sfc, 12, loadLUT());
        sfcR = decodeHuff(AACSeq3(i).chr.sfc, 12, loadLUT());
        S_L = decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook, loadLUT());
        S_R = decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook, loadLUT());
        
        frameF_L = iAACquantizer(S_L, [AACSeq3(i).chl.G ; sfcL'], AACSeq3(i).chl.G, AACSeq3(i).frameType);
        frameF_R = iAACquantizer(S_R, [AACSeq3(i).chr.G ; sfcR'], AACSeq3(i).chr.G, AACSeq3(i).frameType);
        frameF = [iTNS(frameF_L, AACSeq3(i).frameType, AACSeq3(i).chl.TNScoeffs) iTNS(frameF_R, AACSeq3(i).frameType, AACSeq3(i).chr.TNScoeffs)] ;
    end
    
    frameT = iFilterbank( frameF, AACSeq3(i).frameType, AACSeq3(i).winType );
    
    %Add overlapping parts of frames
    channel_1 = [ channel_1 ; prevFrameT(1025:2048,1)+frameT(1:1024,1) ];
    channel_2 = [ channel_2 ; prevFrameT(1025:2048,2)+frameT(1:1024,2) ];
    
    prevFrameT = frameT;

end

%Exclude last frame
%channel_1 = [ channel_1 ; frameT(1025:2048,1) ];
%channel_2 = [ channel_2 ; frameT(1025:2048,2) ];

y = [channel_1 channel_2 ];

%Write decoded audio
audiowrite(fNameOut, y, 48000) ;

if nargout == 1
    x = y;
end

end

