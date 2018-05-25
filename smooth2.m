function x = smooth2 (x, n)
for i = 1:size(x,2)
    g = get_contig_groups(find(~isnan(x(:,i))));
    for j = 1:length(g)
        a = g{j}(1);
        b = g{j}(end);
        x(a:b,i) = smooth(x(a:b,i),n);
    end
end

