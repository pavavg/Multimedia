function SNR = demoAAC1(fNameIn, fNameOut)
AACSeq1 = AACoder1(fNameIn);
y2=iAACoder1(AACSeq1, fNameOut);
y1 = audioread(fNameIn);
%y2 = audioread(fNameOut);

noise1 = y1(1025:size(y2,1)+1024, 1 ) -y2 (:,1);
noise2 = y1(1025:size(y2,1)+1024, 2 ) -y2 (:,2);
SNR =  [ snr(y1(1025:size(y2,1)+1024,1), noise1) snr(y1(1025:size(y2,1)+1024,2), noise2) ]; 
10*log10(var(y1(1:size(y2,1),1))/var(noise1))
plot(noise1)
figure
plot(noise2)
end

