function frameT = iFilterbank(frameF, frameType, winType)
    
    N = 2048 ;
if strcmp(frameType, 'OLS')
    X1 = imdct4 (frameF(:,1) );
    X2 = imdct4 (frameF(:,2) );
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
        
        S1 = X1.*W;
        S2 = X2.*W;
        
    elseif strcmp(winType, 'SIN')
        
        Wleft = zeros(M,1);
        Wright = zeros(M,1);
        
        for n=1:M
            Wleft(n) = sin(pi *( n +0.5) /N ) ;
            Wright(n) = sin(pi *( M+n +0.5) /N ) ;
        end
        W = [Wleft ; Wright];
        
        S1 = X1.*W;
        S2 = X2.*W;
        
    end
    frameT = [S1 S2] ;

elseif strcmp(frameType, 'LSS')
     X1 = imdct4 (frameF(:,1) );
     X2 = imdct4 (frameF(:,2) );
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
        
        S1 = X1.*W;
        S2 = X2.*W;
        
        
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
        
        S1 = X1.*W;
        S2 = X2.*W;
    end
    frameT = [S1 S2] ;
elseif strcmp(frameType, 'LPS')
     X1 = imdct4 (frameF(:,1) );
     X2 = imdct4 (frameF(:,2) );
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
        
        S1 = X1.*W;
        S2 = X2.*W;
        
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
        
        S1 = X1.*W;
        S2 = X2.*W;
        
    end
    frameT = [S1 S2] ;    
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
    
    F1 = zeros(128,8);
    F2 = zeros(128,8);
    S1 = zeros(256, 8) ;
    S2 = zeros(256, 8) ;
    S1_all =[];
    S2_all =[];
    index = 1;
    for i =1:8
        F1(:,i) = frameF(index: index+127 , 1);
        F2(:,i) = frameF(index: index+127 , 2);
        
        S1(:,i) = imdct4( F1(:,i) ) .* W ;
        S2(:,i) = imdct4( F2(:,i) ) .* W ;
        
        S1_all = [S1_all ; S1(:,i)];
        S2_all = [S2_all ; S2(:,i)];
        
        index = index +128 ;
    end
    frameT = [S1_all S2_all] ;
end

end

