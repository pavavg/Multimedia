function frameFout = iTNS(frameFin, frameType, TNScoeffs)

if ~strcmp(frameType, 'ESH')
    enum = [1 TNScoeffs] ;
    denom = 1 ;
    
    frameFout = filter(denom, enum, frameFin) ;  
    
else
    frameFout = zeros(128,8);
    for i = 1:8
        enum = [1 ;TNScoeffs(:,i)] ;
        denom = 1 ;

        frameFout(:,i) = filter(denom, enum, frameFin(:,i)) ;
    
    end
end



end

