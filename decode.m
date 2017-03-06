function cdwrd_est = decode(rcvd_dem,n,k,H,I,type)
    if (strcmp(type,'repetition'))
        if (mean(rcvd_dem) > 0.5)
            cdwrd_est = ones(1,n); 
        else
            cdwrd_est = zeros(1,n);
        end
        
    elseif (strcmp(type,'none'))
        cdwrd_est = rcvd_dem;
        
    elseif(strcmp(type,'Hamming'))
        e = rcvd_dem*H.';
        e1 = mod(e,2);
        
        rcvd_corrected = rcvd_dem;
        
        % we can only correct errors of weight 1
        if (sum(e1) == 1)
            location = ismember(H',e1,'rows');
            rcvd_corrected(location) = mod(rcvd_corrected(location)+1,2);
        end
        cdwrd_est = rcvd_corrected;
    end
end
    