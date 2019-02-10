function frameFout = iTNS(frameFin, frameType, TNScoeffs)

if ~strcmp(frameType, 'ESH')
    enum = TNScoeffs ;
    denom = 1 ;
    
    frameFout = filter(denom, enum, frameFin) ;  
    
else
    for i = 1:8
        enum = TNScoeffs(:,i) ;
        denom = 1 ;

        frameFout(i) = filter(denom, enum, frameFin(:,i)) ;
    
    end
end



end

