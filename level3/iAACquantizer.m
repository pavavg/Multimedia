function frameF = iAACquantizer(S, sfc, G, frameType)

load('TableB219.mat') ;

if ~strcmp (frameType, 'ESH')
    a = zeros(size(B219a , 1), 1) ;
    a(1) = G ;
    Xhat(1) = sign(S(1))*(abs(S(1))).^(4/3) * 2^(a(1)/4) ;

    for i = 2:size(sfc) 
        a(i) = sfc(i) - a(i-1) ;
        Xhat(i) = sign(S(i))*(abs(S(i)))^(4/3) * 2^(a(i)/4) ;
    end
    
    
else
    a = zeros(size(B219b , 1), 8) ; 
    a(1, :) = G ;
    Xhat(1, :) = sign(S(1, :))*(abs(S(1, :)))^(4/3) * 2^(a(1, :)/4) ;
    
    for f = 1:8
        for i = 2:size(sfc(:, f))
            a(i, f) = sfc(i, f) - a(i-1, f) ;
            Xhat(i, f) = sign(S(i, f))*(abs(S(i, f)))^(4/3) * 2^(a(i, f)/4) ;
        end   
    end
    
end

frameF = Xhat ;

end