function pvalue = multiary_task_chance_levels (total_number_of_trials, total_number_of_corrects, chance_probability)
%{
Tells you the probability that by pure chance you can get the total_number_of_corrects or higher
Consequently tells you whether your decision maker is random
%}


n = total_number_of_trials;
k = [ceil(total_number_of_corrects):n];
p = chance_probability;
q = 1 - p;
l = n-k;
logfacn = logfac(n);
lp = log(p);
lq = log(q);

logp = nan(1,length(k));

for i = 1:length(k)
    logp(i) = logfacn - logfac(l(i)) - logfac(k(i)) + k(i)*lp + (l(i))*lq;
end

pvalue = sum(exp(logp));


function lf = logfac (n)
%log factorial
%log(n!) = \sum_{i=1}^{n} log(i)
lf = sum(log(1:n));
