function rcvd_dem = demodulate(rcvd,symb,n) 
    rcvd_dem = zeros(1,2*length(symb));
    rcvd_dem(1:2:end) = (real(rcvd)<0);
    rcvd_dem(2:2:end) = (imag(rcvd)<0);
    rcvd_dem = rcvd_dem(1:n);
end