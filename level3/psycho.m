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
        hannWindow = 0.5 -0.5* cos( pi*(i-1+0.5) / 1024 );
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
    c_w = sqrt( ( r.*cos(f)-rpred.*cos(fpred) ).^2 + ( r.*sin(f)-rpred.*sin(f) ).^2 ) ./ (r + rpred) ;
    
    e = zeros (69,1);
    c = zeros (69,1);
    
    %Query 5
    for i = 1:69
        start = longTable(i,2) +1;
        finish = longTable(i,3) +1;
        e(i) = sum ( r(start:finish).^2 );
        c(i) = sum( c_w(start:finish).* r(start:finish).^2  );
    end
    
    %Query 6
    ecb = zeros(69,1);
    ct = zeros(69,1);
    en = zeros(69,1);
    for b = 1:69
        ecb(b) = sum ( e(1:69) .* SpreadFunLong(1:69, b) );
        ct(b) = sum ( c(1:69) .* SpreadFunLong(1:69, b) );
        en(b) = ecb(b) / sum ( SpreadFunLong(1:69, b) ) ;
    end
    cb = ct ./ ecb ;
    
    %Query 7
    tb = -0.299 - 0.43 * log(cb) ;
    
    %Query 8
    SNR = tb *18 + (1-tb)* 6 ;
    
    %Query 9
    bc = 10 .^ (-SNR/10);
    
    %Query 10
    nb = en .* bc;
    
    %Query 11
    qthr = eps() * 1024 * 10.^(longTable(:,6)/10) ;
    npart = max( nb, qthr );
    
    %Query 12
    SMR = e ./npart ;
    

else  
    
    %Query 1
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
    
    %Query 2
    Sub = zeros(256,8) ;
    Subprev1 = zeros(256,8) ;
    Sw = zeros(256,8) ;
    Swprev1 = zeros(256,8) ;
    Y = zeros(128,8);
    Yprev1 = zeros(128,8);
    r = zeros(128,8);
    f = zeros(128,8);
    rprev1 = zeros(128,8);
    fprev1 = zeros(128,8);
    for i =1:8
        index = 449 + i * 128;
        Sub(:,i) = frameT(index:index+255);
        Subprev1(:,i) = frameTprev1(index:index+255);
        
        for n =1:256
            hannWindow = 0.5 -0.5* cos( pi*(n-1+0.5) / 128 );
            Sw(n,i) = Sub(n,i)*hannWindow ;
            Swprev1(n,i) = Subprev1(n,i)*hannWindow ;
        end
        
        Y(:,i) = fft(Sw(:,i), 128);
        Yprev1(:,i) = fft(Swprev1(:,i), 128);
        
        r(:,i) = abs( Y(:,i) );
        f(:,i) = angle( Y(:,i) );
        
        rprev1(:,i) = abs( Yprev1(:,i) );
        fprev1(:,i) = angle( Yprev1(:,i) );
    end

    %Query 3
    rpred = zeros(128,8);
    fpred = zeros(128,8);
    
    rpred(:,1) = 2*rprev1(:,8); - rprev1(:,7) ;
    fpred(:,1) = 2*fprev1(:,8); - fprev1(:,7) ;
    
    rpred(:,2) = 2*r(:,1); -  rprev1(:,8);
    fpred(:,2) = 2*f(:,1); -  fprev1(:,8);
    
    for i = 3:8
        rpred(:,i) = 2*r(:,i-1) - r(:, i-2);
        fpred(:,i) = 2*f(:,i-1) - f(:, i-2);
    end
    
    %Query 4
    c_w = sqrt( ( r.*cos(f)-rpred.*cos(fpred) ).^2 + ( r.*sin(f)-rpred.*sin(f) ).^2 ) ./ (r + rpred) ;
    
    %Query 5
    e = zeros(42,8);
    c = zeros(42,8);
    for i=1:8
        for b = 1:42
            start = shortTable(b,2) +1;
            finish = shortTable(b,3) +1;
            e(b,i) = sum ( r(start:finish,i).^2 );
            c(b,i) = sum ( c_w(start:finish,i) .* r(start:finish,i).^2 );
        end
    end
    
    %Query 6
    ecb = zeros(42,8);
    ct = zeros(42,8);
    en = zeros(42,8);
    for i = 1:8
        for b = 1:42
            ecb(b,i) = sum ( e(1:42,i) .* SpreadFunShort(1:42, b) );
            ct(b,i) = sum ( c(1:42,i) .* SpreadFunShort(1:42, b) );
            en(b,i) = ecb(b,i) / sum ( SpreadFunShort(1:42, b) ) ;
        end
    end
    cb = ct ./ ecb ;
    
    %Query 7
    tb = -0.299 - 0.43*log(cb) ;
    
    %Query 8
    SNR = tb*18 + (1-tb)*6;
    
    %Query 9
    bc = 10 .^(-SNR/10) ;
    
    %Query 10
    nb = en .* bc;
    
    %Query 11
    qthr = zeros(42,8);
    for i = 1:8
        qthr(:,i) = eps() * 128 * 10.^(shortTable(:,6)/10) ;
    end
    npart = max(nb , qthr);
    
    %Query 12
    SMR = e ./npart;
    
end
    
end

