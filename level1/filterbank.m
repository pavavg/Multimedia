function frameF = filterbank(frameT, frameType, winType)
N = 2048 ; %frame Size

% Create window W and calculate the new signal
if strcmp(frameType, 'OLS')
    
    M = N/2;
    if strcmp(winType, 'KBD' )
             
        kais = kaiser(M+1,6*pi);
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
            Wleft(n) = sin(pi *( n-1 +0.5) /N ) ;
            Wright(n) = sin(pi *( M+n-1 +0.5) /N ) ;
        end
        W = [Wleft ; Wright];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
    end

elseif strcmp(frameType, 'LSS')
    
    if strcmp(winType, 'KBD' )
        M = N/2;
        kaisLeft = kaiser(M+1,6*pi);
        Wleft = zeros(M,1);
        kaisLeft_sum = sum(kaisLeft);
        for n=1:M
            Wleft(n) = sqrt( sum(kaisLeft(1:n)) / kaisLeft_sum );
        end
        
        Wright1 = ones(448,1);
        
        kaisRight2= kaiser(128+1,4*pi);
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
            Wleft(n) = sin(pi *( n-1 +0.5) /N ) ;
        end
        
        Wright1 = ones(448,1);
        
        Wright2 = zeros(128,1);
        
        for n=1:128
            Wright2(n) = sin(pi *(128+ n-1 +0.5) /256 ) ;
        end
        
        Wright3 = zeros(448,1);
        
        W = [Wleft ; Wright1 ; Wright2 ; Wright3 ];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
    end
    
elseif strcmp(frameType, 'LPS')
    
    if strcmp(winType, 'KBD' )
        Wleft1 = zeros(448,1);
        
        kaisleft2= kaiser(128+1,4*pi);
        Wleft2 = zeros(128,1);
        kaisleft2_sum = sum(kaisleft2);
        for n=1:128
            Wleft2(n) = sqrt( sum(kaisleft2(1:n)) / kaisleft2_sum );
        end
        
        Wleft3 = ones(448,1);
        
        M = N/2;
        kaisRight = kaiser(M+1,6*pi);
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
            Wleft2(n) = sin(pi *( n-1+0.5) /256 ) ;
        end
        
        Wleft3 = ones(448,1);
        
        M = N/2;
        Wright = zeros(M,1);
        
        for n=1:M
            Wright(n) = sin(pi *(M+ n-1 +0.5) /N ) ;
        end
        
        W = [Wleft1 ; Wleft2 ; Wleft3 ;Wright  ];
        
        Snew1 = W.* frameT(:,1) ;
        Snew2 = W.* frameT(:,2) ;
        
    end
        
else
    %Create W for subframes
    if strcmp(winType, 'KBD' )
        kais = kaiser(128+1,4*pi);
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
            Wleft(n) = sin(pi *( n-1 +0.5) /256 ) ;
            Wright(n) = sin(pi *( 128+n-1 +0.5) /256 ) ;
        end
        W = [Wleft ; Wright];
    end
    %Apply W to each subframe
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

%Calculate MDCT
if strcmp(frameType, 'LPS') || strcmp(frameType, 'OLS') || strcmp(frameType, 'LSS')
    
    X1 = mdct4(Snew1);
    X2 = mdct4(Snew2);
    
    frameF = [X1 X2] ;

else
    
    frameF1 =[] ;
    frameF2 =[] ;
    for i = 1:8
        X1 = mdct4(Snew1(:,i));
        X2 = mdct4(Snew2(:,i));
        frameF1 = [frameF1 ; X1];
        frameF2 = [frameF2 ; X2];
        
    end
        
        frameF = [frameF1 frameF2] ;
end
    

end