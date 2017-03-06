function symb = modulate(cdwrd,Es)
    if (mod(length(cdwrd),2) == 1)
        cdwrd = [cdwrd 0];
    end
    symb = sqrt(Es/2)*(1-2*cdwrd(1:2:end))+1i*sqrt(Es/2)*(1-2*cdwrd(2:2:end)); 
end