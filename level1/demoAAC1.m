function SNR = demoAAC1(fNameIn, fNameOut)
AACSeq1 = AACoder1(fNameIn);
iAACoder1(AACSeq1, fNameOut);
y1 = audioread(fNameIn);
y2 = audioread(fNameOut);

noise1 = y1(1:size(y2,1), 1 ) -y2 (:,1);
noise2 = y1(1:size(y2,1), 2 ) -y2 (:,2);
SNR =  [ snr(y1(1:size(y2,1),1), noise1) snr(y1(1:size(y2,1),2), noise2) ]; 

plot(noise1)
end

