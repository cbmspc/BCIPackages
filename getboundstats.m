function boundstats = getboundstats(x, Fs, bounds_sec)
for i = size(bounds_sec,1):-1:1
    a = round(bounds_sec(i,1)*Fs+1);
    b = round(bounds_sec(i,2)*Fs);
    boundstats(i).t_range = bounds_sec(i,1:2);
    boundstats(i).s_range = [a b];
    boundstats(i).min = min(x(a:b,:));
    boundstats(i).max = max(x(a:b,:));
    boundstats(i).mean = mean(x(a:b,:));
    boundstats(i).median = median(x(a:b,:));
    boundstats(i).std = std(x(a:b,:));
    boundstats(i).robuststd = robuststd(x(a:b,:));
end
