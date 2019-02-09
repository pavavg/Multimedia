function frameF = filterbank(frameT, frameType, winType)
N = 2048 ;
if strcmp(frameType, 'OLS')
    
    M = N/2;
    if strcmp(winType, 'KBD' )
             
        kais = kaiser(M,6);
        Wleft = zeros(M,1);
        Wright = zeros(M,1);
        kais_sum = sum(kais);
        for n=1:M
            Wleft(n) = sqrt( sum(kais(1:n)) / kais_sum );
            Wright(n) = sqrt( sum(kais(1: M+1-n )) / kais_sum);
        end
        W = [Wleft ; Wright];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
        
    elseif strcmp(winType, 'SIN')
        
        Wleft = zeros(M,1);
        Wright = zeros(M,1);
        
        for n=1:M
            Wleft(n) = sin(pi *( n +0.5) /N ) ;
            Wright(n) = sin(pi *( M+n +0.5) /N ) ;
        end
        W = [Wleft ; Wright];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
    end

elseif strcmp(frameType, 'LSS')
    
    if strcmp(winType, 'KBD' )
        M = N/2;
        kaisLeft = kaiser(M,6);
        Wleft = zeros(M,1);
        kaisLeft_sum = sum(kaisLeft);
        for n=1:M
            Wleft(n) = sqrt( sum(kaisLeft(1:n)) / kaisLeft_sum );
        end
        
        Wright1 = ones(448,1);
        
        kaisRight2= kaiser(128,6);
        Wright2 = zeros(128,1);
        kaisRight2_sum = sum(kaisRight2);
        for n=1:128
            Wright2(n) = sqrt( sum(kaisRight2(1:129-n)) / kaisRight2_sum );
        end
        
        Wright3 = zeros(448,1);
        
        W = [Wleft ; Wright1 ; Wright2 ; Wright3 ];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
        
    elseif strcmp(winType, 'SIN')
        M = N/2;
        Wleft = zeros(M,1);
        for n=1:M
            Wleft(n) = sin(pi *( n +0.5) /N ) ;
        end
        
        Wright1 = ones(448,1);
        
        Wright2 = zeros(128,1);
        
        for n=1:128
            Wright2(n) = sin(pi *(128+ n +0.5) /256 ) ;
        end
        
        Wright3 = zeros(448,1);
        
        W = [Wleft ; Wright1 ; Wright2 ; Wright3 ];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
    end
    
elseif strcmp(frameType, 'LPS')
    
    if strcmp(winType, 'KBD' )
        Wleft1 = zeros(448,1);
        
        kaisleft2= kaiser(128,6);
        Wleft2 = zeros(128,1);
        kaisleft2_sum = sum(kaisleft2);
        for n=1:128
            Wleft2(n) = sqrt( sum(kaisleft2(1:n)) / kaisleft2_sum );
        end
        
        Wleft3 = ones(448,1);
        
        M = N/2;
        kaisRight = kaiser(M,6);
        Wright = zeros(M,1);
        kaisRight_sum = sum(kaisRight);
        for n=1:M
            Wright(n) = sqrt( sum(kaisRight(1:M+1-n)) / kaisRight_sum );
        end
        
        W = [Wleft1 ; Wleft2 ; Wleft3 ;Wright  ];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
    elseif strcmp(winType, 'SIN')

        Wleft1 = zeros(448,1);

        Wleft2 = zeros(128,1);

        for n=1:128
            Wleft2(n) = sin(pi *( 128 +0.5) /256 ) ;
        end
        
        Wleft3 = ones(448,1);
        
        M = N/2;
        Wright = zeros(M,1);
        
        for n=1:M
            Wright(n) = sin(pi *(M+ n +0.5) /N ) ;
        end
        
        W = [Wleft1 ; Wleft2 ; Wleft3 ;Wright  ];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
    end
        
else
    
    if strcmp(winType, 'KBD' )
        kais = kaiser(128,6);
        Wleft = zeros(128,1);
        Wright = zeros(128,1);
        kaisLeft_sum = sum(kais);
        for n=1:128
            Wleft(n) = sqrt( sum(kais(1:n)) / kaisLeft_sum );
            Wright(n) = sqrt( sum(kais(1: 128+1-n )) / kaisLeft_sum);
        end
        W = [Wleft ; Wright];
        
    elseif strcmp(winType, 'SIN' )
        Wleft = zeros(128,1);
        Wright = zeros(128,1);
        
        for n=1:128
            Wleft(n) = sin(pi *( n +0.5) /256 ) ;
            Wright(n) = sin(pi *( 128+n +0.5) /256 ) ;
        end
        W = [Wleft ; Wright];
    end
    
    S1 = zeros(256,8) ;
    S2 = zeros(256,8) ;
    Snew1 = zeros(256,8) ;
    Snew2 = zeros(256,8) ;
    for i =1:8
        index = 449 + (i-1) * 128;
        S1(:,i) = frameT(index:index+255, 1);
        S2(:,i) = frameT(index:index+255, 2);
        Snew1(:,i) = W .* S1(:,i) ;
        Snew2(:,i) = W .* S2(:,i) ;
    end
end


if strcmp(frameType, 'LPS') || strcmp(frameType, 'OLS') || strcmp(frameType, 'LSS')
    X1 = zeros(1024,1);
    X2 = zeros(1024,1);
    n0 = (1024+1) /2 ;
    for k = 1: 1024
        for n = 1:2048
            X1(k) = X1(k) + Snew1(n) * cos( (2*pi/2048) * (n+n0) * (k+0.5) );
            X2(k) = X2(k) + Snew2(n) * cos( (2*pi/2048) * (n+n0) * (k+0.5) );
        end
    end
    
    frameF = 2 * [X1 X2] ;

else
    
    frameF1 =[] ;
    frameF2 =[] ;
    X1 = zeros(128,1);
    X2 = zeros(128,1);
    n0 = (128+1) /2 ;
    for k = 1: 128
        for n = 1:256
            X1(k) = X1(k) + Snew1(n,i) * cos( (2*pi/256) * (n+n0) * (k+0.5) );
            X2(k) = X2(k) + Snew2(n,i) * cos( (2*pi/256) * (n+n0) * (k+0.5) );
        end
    end
    frameF = [X1 X2];
        
    for i = 2:8
        X1 = zeros(128,1);
        X2 = zeros(128,1);
        n0 = (128+1) /2 ;
        for k = 1: 128
            for n = 1:256
                X1(k) = X1(k) + Snew1(n,i) * cos( (2*pi/256) * (n+n0) * (k+0.5) );
                X2(k) = X2(k) + Snew2(n,i) * cos( (2*pi/256) * (n+n0) * (k+0.5) );
            end
        end
        
        frameF = cat(3,frameF,[X1 X2]);
        
    end
        
        frameF = 2 * frameF ;
end
    
end