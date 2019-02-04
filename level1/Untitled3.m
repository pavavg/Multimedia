N= 2048;
M = N/2;
 kais = kaiser(128+1,4*pi);
        Wleft = zeros(128,1);
        Wright = zeros(128,1);
        kaisLeft_sum = sum(kais);
        for n=1:128
            Wleft(n) = sqrt( sum(kais(1:n)) / kaisLeft_sum );
            Wright(n) = sqrt( sum(kais(1: 128+1-n )) / kaisLeft_sum);
        end
        W = [Wleft ; Wright];
plot(W)