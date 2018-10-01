function [sol,a,b] = scaledPartialGauss(n,A0,b0)
%SCALEDPARTIALGAUSS Gaussian Elimination with scaled partial pivoting
%   Takes size of matrix (nxn), coefficient matrix, b vector.
%   Returns solutions and edited matric and b vector.

s = zeros(n,1); %allocate s
l = (linspace(1,n,n))';
a = A0;
b= b0;

%Foreward Elimination
for i = 1:n
    l(i) = i;
    smax = 0;
    for j = 1:n
        smax = max([smax,abs(a(i,j))]);
    end
    s(i) = smax;
end
for k = 1:(n-1)
    rmax = 0;
    for i = k:n
        r = abs(a(l(i),k)/s(l(i)));
        if r > rmax
            rmax = r;
            j = i;
        end
    end
    hold = l(j);
    l(j) = l(k);
    l(k) = hold;
    for i = (k+1):n
        xmult = a(l(i),k)/a(l(k),k);
        a(l(i),k) = xmult;
        for j = (k+1):n
            a(l(i),j) = a(l(i),j) - (xmult * a(l(k),j));
        end
    end
end

%Foreward Elimination for b
for k = 1:(n-1)
    for  i = (k+1):n
        b(l(i)) = b(l(i)) - (a(l(i),k) * b(l(k)));
    end
end

%Back substitution
x(n) = b(l(n))/a(l(n),n);
for i = (n-1):-1:1
    sum = b(l(i));
    for j = (i+1):n
        sum = sum - (a(l(i),j) * x(j));
    end
    x(i) = sum/a(l(i),i);
end

sol = x;
end

