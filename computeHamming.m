function [G,H,I] = computeHamming(n,k)
    h1 = 1:n;
    for i = 1:n
        if mod(log2(h1(i)),1)==0
            v = h1(log2(i)+1);
            h1(log2(i)+1) =  h1(i);
            h1(i) = v;
        end
    end
 
    H1 = zeros(n-k,n);
    for i=1:n
        H1(:,i) =de2bi(h1(i),n-k)';
    end
   
    A = H1(:,n-k+1:n);
    G1 = [A' eye(k)];
    
    % we Have H and G in systematic form, we need to put H in its
    % 1,2,3,4,5,..... form and G must follow the same column swaps
    [H2,I] = sort(h1,2);
    H = H1(:,I);
    G = G1(:,I);
end
