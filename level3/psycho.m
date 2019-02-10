function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)

    
load('TableB219.mat')
longTable = B219a;
shortTable = B219b;
if ~strcmp (frameType, 'ESH')
    
    %Query 1
    SpreadFunLong = zeros(size(longTable,1));
    for i=1:size(longTable,1)
        for j=1:size(longTable,1)

            if i>=j
                tmpx = 3*(longTable(j,5) -longTable(i,5) );
            else
                tmpx = 1.5*(longTable(j,5) -longTable(i,5) );
            end

            tmpz = 8*min((tmpx-0.5)^2 - 2*(tmpx-0.5),0);
            tmpy = 15.811389 +7.5*(tmpx+0.474) - 17.5* sqrt(1 + (tmpx + 0.474)^2);

            if tmpy< -100
                SpreadFunLong(i,j) = 0;
            else
                SpreadFunLong(i,j) = 10 ^ ( (tmpz+tmpy)/10 );
            end
        end
    end
    
    %Query 2
    Sw = zeros(2048,1) ;
    Swprev1 = zeros(2048,1) ;
    Swprev2 = zeros(2048,1) ;
    for i = 1: 2048
        hannWindow = 0.5 -0.5* cos( pi*(n-1+0.5) / 2048 );
        Sw(i) = frameT(i)*hannWindow ;
        Swprev1(i) = frameTprev1(i)*hannWindow ;
        Swprev2(i) = frameTprev2(i)*hannWindow ;
    end
    
    Y1 = fft(Sw, 1024) ;
    Y2 = fft(Swprev1, 1024) ;
    Y3 = fft(Swprev2, 1024) ;
    
    r = abs(Y1);
    f = angle(Y1);
    
    rprev1 = abs(Y2);
    fprev1 = angle(Y2);
    
    rprev2 = abs(Y3);
    fprev2 = angle(Y3);
    
    %Query 3
    rpred = 2*rprev1 - rprev2;
    fpred = 2*fprev1 - fprev2;
    
    %Query 4
    c = sqrt( ( r.*cos(f)-rpred.*cos(fpred) )^2 + ( r.*sin(f)-rpred.*sin(f) )^2 ) / (r + rpred) ;
    
    
SpreadFunShort = zeros(size(shortTable,1));
for i=1:size(shortTable,1)
    for j=1:size(shortTable,1)

        if i>=j
            tmpx = 3*(shortTable(j,5) -shortTable(i,5) );
        else
            tmpx = 1.5*(shortTable(j,5) -shortTable(i,5) );
        end

        tmpz = 8*min((tmpx-0.5)^2 - 2*(tmpx-0.5),0);
        tmpy = 15.811389 +7.5*(tmpx+0.474) - 17.5* sqrt(1 + (tmpx + 0.474)^2);

        if tmpy< -100
            SpreadFunShort(i,j) = 0;
        else
            SpreadFunShort(i,j) = 10 ^ ( (tmpz+tmpy)/10 );
        end
    end
end


end

