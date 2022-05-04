function ildDB = getILD(irLeft, irRight, Fs)

    %% ERB stuff
    erb = 1:40;
    Q = 9.265;
    L = 24.7;
    gamma = 1;
    for i = erb
        fc(i) = Q * L * (exp(i*gamma/Q) - 1);
        fbw(i) = gamma * L * exp(i*gamma/Q);
    end

    %% get magnitude vectors and ERB
    [f,magL] = getMagnitude(irLeft,Fs,'lin');
    [f,magR] = getMagnitude(irRight,Fs,'lin');

    for i = erb
       erb_pwrL(i) = sum(magL(f >= fc(i)-fbw(i)/2 & f < fc(i)+fbw(i)/2).^2);
       erb_pwrR(i) = sum(magR(f >= fc(i)-fbw(i)/2 & f < fc(i)+fbw(i)/2).^2);
    end
    
    idx = find(fc >= 1500);
    ildDB = 10*log10(sum(erb_pwrL(idx))/sum(erb_pwrR(idx)));

end

