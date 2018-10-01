function [x] = naiveGauss(n,a,b,x)
%NAIVEGAUSS Performs Naive Gaussian elimination
%   Takes in matrix size, matrix, results, and empty solutions array.
%   Outputs System solutions as an n element array.

for k = 1:n-1 %for each iteration
    for ind = k+1:n %for each row
        xMult = a(ind,k)./a(k,k); %set multiplier
%         a(ind,k) = xMult; %replace entry with multiplier
        for j = k:n %for each column
            a(ind,j) = a(ind,j) - (xMult .* a(k,j)); %multiply and subtract pivot 
        end
        b(ind) = b(ind) - (xMult .* b(k)); %multiply and subtract pivot
    end
end
x(n) = b(n)./a(n,n); %set nth x
for indi = 1:n-1 %for each row
    indix = n - indi;
    sum = b(indix); %set sum of row
    for j = indix+1:n %for each column
        sum = sum - (a(indix,j) .* x(j)); %sum the known elements
    end
    x(indix) = sum ./ a(indix,indix); %assign the solution for x
end

end

