function frameF = iAACquantizer(S, sfc, G, frameType)

load('TableB219.mat') ;

if ~strcmp (frameType, 'ESH')
    a = zeros(size(B219a , 1), 1) ;
    a(1) = G ;
    start = B219a(1,2)+1;
    finish = B219a(1,3)+1;
    Xhat(start:finish) = sign( S(start:finish) ).*( abs(S(start:finish)) ).^(4/3) * 2^(a(1)/4) ;

    for i = 2:69
        
        start = B219a(i,2)+1;
        finish = B219a(i,3)+1;
        a(i) = sfc(i) + a(i-1) ;
        Xhat(start:finish) = sign(S(start:finish)).* abs( S(start:finish) ).^(4/3) * 2^(a(i)/4) ;
    end
    Xhat = Xhat';
    
else
    a = zeros(size(B219b , 1), 8) ; 
    a(1, :) = G ;
    start = B219b(1,2)+1;
    finish = B219b(1,3)+1;
    Xhat(start:finish, :) = sign(S(start:finish, :)).*(abs(S(start:finish, :))).^(4/3) .* 2.^(a(1, :)/4) ;
    
    for f = 1:8
        for i = 2:42
            start = B219b(i,2)+1;
            finish = B219b(i,3)+1;
            a(i, f) = sfc(i, f) + a(i-1, f) ;
            Xhat(start:finish, f) = sign(S(start:finish, f)).*(abs(S(start:finish, f))).^(4/3) .* 2.^(a(i, f)/4) ;
        end   
    end
    
end

frameF = Xhat ;

end