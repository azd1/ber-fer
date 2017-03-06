function [fer,ber] = computeBFER(n, k, Es, Eb, sigma2, iter, type, file, aff)
    dim = length(sigma2);
    
    fer = zeros(1,dim);
    ber = zeros(1,dim);
    
    H = 0; I = 0;
    try
        load(file);
    catch
        disp 'no file saved before';
        frame_err = zeros(1,dim); % initialize the number of frame errors
        bit_err = zeros(1,dim); % initialize the number of bit errors 
        frame_nb = 0; % initialize the number of codeword transmissions
    end
    frame_nb_init = frame_nb; % offset
    
    if (strcmp(type,'repetition'))
        G = ones(1,n); % generator matrix for repetition code
        H = 0;
        I = 0;
    elseif (strcmp(type,'Hamming'))
        [G,H,I] = computeHamming(n,k);
    elseif (strcmp(type,'Golay'))
        [G,H] = golayMatrices();
        SDT = syndtable(H);
        [k,n] = size(G);
        cdwrds = zeros(1,n);
        msg = dec2bin(0:2^k-1,k) == '1'; % create matrix of all possible messages
        
        for i=1:2^k
        if (~ismember(mod(msg(i,:)*G,2),cdwrds,'rows'))
            cdwrds = [cdwrds; mod(msg(i,:)*G,2)];
        end
        end
        SA = cdwrds;
        eH = zeros(4096,12);
        for i = 1:length( SDT(:,1) )
            eH(i,:) = mod(SDT(i,:)*H.',2);
        end
    end
    
    if (aff)
        PowEff = ones(1,length(sigma2)) * Eb ./ (2*sigma2); % absissa for graphic
        figure
        h = semilogy(10*log10(PowEff),ber);
        grid on
    end

    while(frame_nb < frame_nb_init + iter) % repetition for precision

        frame_nb = frame_nb + 1; % update count of frame simulations
        msg = (rand(1,k)>0.5); % generates k bits uniformly at random
        
        if(~strcmp(type,'none'))
            cdwrd = mod(encode(msg,G),2);
        else 
            cdwrd = msg;
        end
        symb = modulate(cdwrd, Es);
        
        for m = 1:dim    % dim = number of points      
            rcvd = symb + sqrt(sigma2(m))*(randn(size(symb))+1i*randn(size(symb))); % channel noise
            rcvd_dem = demodulate(rcvd,symb,n);
            if strcmp(type,'Golay')
                cdwrd_est = decodeGolay(rcvd_dem,H,SDT,eH);
            else 
                cdwrd_est = decode(rcvd_dem,n,k,H,I,type);
            end
            
            bit_err(m) = bit_err(m) + (sum(cdwrd~=cdwrd_est)); % updates the bit errors
            frame_err(m) = frame_err(m) + (sum(cdwrd~=cdwrd_est)>0); % updates the frame errors
        end
        fer = frame_err / frame_nb; % updates frame error rate
        ber = bit_err / (frame_nb * n); % updates bit error rate
        if (mod(frame_nb,100) == 0)
             save(file,'bit_err','frame_err','frame_nb');
             if(aff)
                 set(h,'XData',10*log10(PowEff),'YData',ber);
                 drawnow
             end
        end
    end
    disp(['Done! Number of frames computed: ', num2str(frame_nb)]);

end

