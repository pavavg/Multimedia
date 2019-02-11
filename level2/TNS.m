function [frameFout, TNScoeffs] = TNS(frameFin, frameType)

% QUERY 1

% read bands
load('TableB219.mat') ;
smallBand = B219b ;
bigBand = B219a;
    
% normalization coefficient Sw -> band
if ~strcmp (frameType, 'ESH')
    Sw = zeros(1024,1) ;
    Xw = zeros(1024,1) ;
    
    for i = 1:69
        startB = bigBand(i,2) + 1 ;
        endB = bigBand(i,3) + 1 ;
        
        elems = frameFin(startB:endB) ;
        % energy of a band
        P = sum(elems.^2) ;
    
        Sw(startB:endB) = sqrt(P) ;
    end
    
    % smoothen
    for k = 1023:-1:1 
        Sw(k) = (Sw(k) + Sw(k+1))/2 ;
    end
    
    for k = 2:+1:1024
        Sw(k) = (Sw(k) - Sw(k-1))/2 ;
    end
    
    % normalize 
    Xw = frameFin./Sw ;     
    
else 
    Sw = zeros(128,8) ;
    Xw = zeros(128,8);
    
    % subframes
    for i = 1:8
        for j = 1:42 
        startB = smallBand(j,2) + 1 ;
        endB = smallBand(j,3) + 1 ;
        
        elems = frameFin(startB:endB, i) ;
    
        % energy of a band
        P = sum(elems.^2) ;
    
        Sw(startB:endB, i) = sqrt(P) ;
        end
        
        % smoothen
        for k = 127:-1:1 
            Sw(k,i) = (Sw(k,i) + Sw(k+1,i))/2 ;
        end

        for k = 2:+1:128
            Sw(k,i) = (Sw(k,i) - Sw(k-1,i))/2 ;
        end
        
        % normalize
        Xw(:,i) = frameFin(:,i)./Sw(:,i) ;
    
    end
end

% QUERY 2
   
if ~strcmp (frameType, 'ESH')
    a = lpc(Xw, 4) ;    
    
else
    a = zeros(5,8);
    for i = 1:8
        
        a(:,i) = lpc(Xw(:,i), 4);
    end
    
end

% QUERY 3

% aQ = a quantisized
if ~strcmp (frameType, 'ESH')
    
    
    % N = 2^4 =16
    aQ = (floor(a(2:5)*10))/10 + 0.05 ;
    aQ = min(aQ, +0.75) ;
    aQ = max(aQ, -0.75) ;    

else
    aQ = zeros(4,8) ;
    
    for i = 1:8          
           aQ(:,i) = (floor(a(2:5,i)*10))/10 + 0.05 ;
           aQ(:,i) = min(aQ(:,i), +0.75) ;
           aQ(:,i) = max(aQ(:,i), -0.75) ;  
    end
    
end

% QUERY 4

if ~strcmp (frameType, 'ESH')
    enum = [1 aQ] ;
    denom = 1 ;
    
    H = filter(enum, denom, frameFin) ;
    % stability check
    % if stableF = 1 ok
    stableF = isstable(denom, enum) ;
    
else
    for i = 1:8
        enum = [1 ;aQ(:,i)] ;
        denom = 1 ;

        H(:,i) = filter(enum, denom, frameFin(:,i)) ;
        % stability check
        % if stableF = 1 ok
        stableF = isstable(denom, enum) ;
    
    end
end

frameFout = H ;
TNScoeffs = aQ ;
    
    
end


 


    
    