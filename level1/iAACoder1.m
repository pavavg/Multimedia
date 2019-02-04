function x = iAACoder1(AACSeq1, fNameOut)

frameNumber = size(AACSeq1 , 2);

channel_1 = [];
channel_2 = [];

frameF = [AACSeq1(1).chl.frameF  AACSeq1(1).chr.frameF];
frameT = iFilterbank( frameF, AACSeq1(1).frameType, AACSeq1(1).winType );
prevFrameT = frameT;

%channel_1 = [ channel_1 ; frameT(1:1024,1)];
%channel_2 = [ channel_2 ; frameT(1:1024,2)];

for i=2:frameNumber
    frameF = [AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF];
    frameT = iFilterbank( frameF, AACSeq1(i).frameType, AACSeq1(i).winType );
    
    
    channel_1 = [ channel_1 ; prevFrameT(1025:2048,1)+frameT(1:1024,1) ];
    channel_2 = [ channel_2 ; prevFrameT(1025:2048,2)+frameT(1:1024,2) ];
    
    prevFrameT = frameT;

end

%channel_1 = [ channel_1 ; frameT(1025:2048,1) ];
%channel_2 = [ channel_2 ; frameT(1025:2048,2) ];

y = [channel_1 channel_2 ];

audiowrite(fNameOut, y, 48000) ;

if nargout == 1
    x = y;
end

end
