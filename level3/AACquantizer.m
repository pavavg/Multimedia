function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)

% QUERY 1
X = frameF ;
MagicNumber = 0.4054 ;
MQ = 8191 ;
load(TableB219.mat) ;


if ~strcmp (frameType, 'ESH')
    a = zeros(size(B219a , 1)) ;    
    a(:) = (16/3)*(log2(max(X).^(3/4))/MQ) ;

    % QUERY 2
    S = sign(X)*round((abs(X)*2.^(-a/4))^(3/4) + MagicNumber) ;
    Xhat = sign(S)*(abs(S))^(4/3) * 2.^(a/4) ;

    noBands = size(B219a, 1) ;
    P = zeros(noBands,1);
    Pe = zeros(noBands,1);
    for i = 1:noBands
        P(i)  = sum(X(B219a(i, 2): B219a(i, 3)).^2) ;
        Pe(i) = sum(X((B219a(i, 2): B219a(i, 3)) - Xhat(B219a(i, 2): B219a(i, 3)) ).^2) ;
    end

    T = P./SMR ; 

    for i = 1:noBands
        while  Pe(i) < T(i) 
            a(i) = a(i) + 1 ;
            % i+1 + i ????
            if  i <noBands && max(abs( a(i+1) - a(i) )) > 60 
                 break;
            end    
            
            for j = B219a(i, 2) : B219a(i, 3)
                S(j) = sign(X(j))*round((abs(X(j))*2^(-a(i)/4))^(3/4) + MagicNumber) ;
                Xhat = sign(S(j))*(abs(S(j)))^(4/3) * 2^(a(i)/4) ;
            end
            
            Pe(i) = sum(X((B219a(i, 2): B219a(i, 3)) - Xhat(B219a(i, 2): B219a(i, 3)) ).^2) ;
        end        
    end
    
    G = a(1) ;    
    sfc = a(2:end) - a(1:end-1) ;
    sfc = [G ;sfc] ;
else    
    a = zeros(size(B219b , 1), 8) ;    
    for f = 1:8             
        a(:, f) = (16/3)*(log2(max(X(:, f)).^(3/4))/MQ) ;
    end


    % QUERY 2
    S = sign(X)*round((abs(X)*2.^(-a/4))^(3/4) + MagicNumber) ;
    Xhat = sign(S)*(abs(S))^(4/3) * 2.^(a/4) ;
    
    noBands = size(B219b, 1) ;
    % subframes
    Pe = zeros(noBands,8);
    P = zeros(noBands,8);
    T = zeros(noBands,8);
    for f = 1:8
    
        
        
        for i = 1:noBands
            P(i, f)  = sum(X(B219a(i, 2): B219a(i, 3), f).^2) ;
            Pe(i, f) = sum(X(B219a(i, 2): B219a(i, 3), f) - Xhat(B219a(i, 2): B219a(i, 3), f ).^2) ;
        end
        
       T(:, f) = P(:, f)./SMR(:, f) ; 
        
       
    for i = 1:noBands
        while  Pe(i, f) < T(i, f) 
            a(i, f) = a(i, f) + 1 ;
            % i+1 + i ????
            if  i <noBands && max(abs( a(i+1, f) - a(i, f) )) > 60 
                 break;
            end    
            
            for j = B219b(i, 2) : B219b(i, 3)
                S(j, f) = sign(X(j, f))*round((abs(X(j, f))*2^(-a(i, f)/4))^(3/4) + MagicNumber) ;
                Xhat(j, f) = sign(S(j, f))*(abs(S(j, f)))^(4/3) * 2^(a(i, f)/4) ;
            end
            
            Pe(i, f) = sum(X(B219b(i, 2): B219b(i, 3), f) - Xhat(B219b(i, 2): B219b(i, 3), f).^2) ;
        end   
        
        
    end
       
       
    end    
    G = a(1,:) ;
    sfc = a(2:end, :) - a(1:end-1, :) ;
    sfc = [G ;sfc] ;
    
end

        
end
