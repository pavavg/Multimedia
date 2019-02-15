function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
tic
AACSeq3 = AACoder3(fNameIn, frameAACoded);    %Read frame Structs
toc

tic
y2 = iAACoder3(AACSeq3, fNameOut);    %Decode
toc

y1 = audioread(fNameIn);    %Read original audio  

noise1 = y1(1025:size(y2,1)+1024, 1 ) -y2 (:,1);
noise2 = y1(1025:size(y2,1)+1024, 2 ) -y2 (:,2);

%Calculate SNR for each channel
SNR =  [ snr(y1(1025:size(y2,1)+1024,1), noise1) snr(y1(1025:size(y2,1)+1024,2), noise2) ]; 

info1 = audioinfo(fNameIn) ;
info2 = audioinfo(fNameOut) ;

bitrate1 = info1.BitsPerSample * info1.SampleRate * 2 ;
bitrate2 = info2.BitsPerSample * info2.SampleRate * 2 ;

bitrate = bitrate2 ;
compression = bitrate1 / bitrate2 ;

end

