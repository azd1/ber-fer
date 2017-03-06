function cdwrd_est = decodeGolay(rcvd_dem,H,SDT,eH)
    
    cdwrd_est=rcvd_dem;
    if (~mod(rcvd_dem *H.',2)==0)
        for i=1:4096
            if sum(mod(rcvd_dem *H.',2) == eH(i,:))== 12
               cdwrd_est = mod(rcvd_dem+SDT(i,:),2);
               break
            end
        end
    end    

end
