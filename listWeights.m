function seq = listWeights(n,w)
    if (n == w && n > 0) % case there is no left space for more zeros
        seq = ones(1,n);
    elseif (n == 1) % case 1 elements to add, depends on the left weight
        if (w == 0)
            seq = 0;
        else
            seq = 1;
        end
    elseif (w > 0) % case n elements to add with positive weight
        [a,~] = size(listWeights(n-1,w));
        [c,~] = size(listWeights(n-1,w-1));
        seq = [zeros(a,1) listWeights(n-1,w)
               ones(c,1) listWeights(n-1,w-1)];
    else % case n elements to add with zero weight
        seq = zeros(1,n);
    end
end