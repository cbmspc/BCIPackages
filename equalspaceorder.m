function x2 = equalspaceorder (x)
N = size(x,1);
x2 = x;
k1 = 0;
k2 = round(N/2);
for i = 1:N
    if mod(i,2) == 1
        k1 = k1 + 1;
        x2(i,:) = x(k1,:);
    else
        k2 = k2 + 1;
        x2(i,:) = x(k2,:);
    end
end
