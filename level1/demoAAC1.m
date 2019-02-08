function SNR = demoAAC1(fNameIn, fNameOut)
AACSeq1 = AACoder1(fNameIn);    %Read frame Structs
y2=iAACoder1(AACSeq1, fNameOut);    %Decode
y1 = audioread(fNameIn);    %Read original audio

%Calculate noise for each channel
noise1 = y1(1025:size(y2,1)+1024, 1 ) -y2 (:,1);
noise2 = y1(1025:size(y2,1)+1024, 2 ) -y2 (:,2);

%Calculate SNR for each channel
SNR =  [ snr(y1(1025:size(y2,1)+1024,1), noise1) snr(y1(1025:size(y2,1)+1024,2), noise2) ]; 

end

