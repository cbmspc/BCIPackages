% Unweighted K-nearest neighbor probability
% Returns the probability that "wantlabel" is the predicted label of "test"
% data given "train" training data and "K" neighbors

function Probability = dataproc_func_knnprob (test, train, label, K, wantlabel)

labelnames = unique(label);
Nclass = length(labelnames);
Ntest = size(test,1);
Ntrain = size(train,1);
Fdim = size(train,2);

wantid = find(labelnames == wantlabel);

if Fdim ~= size(test,2)
    error('Feature dimension mismatch.');
end

if length(who('K')) == 0 || length(K) == 0 || K == 0
    K = ceil(Ntrain^(0.25));
elseif K < 0
    K = -K;
    K = ceil(K);
else
    K = ceil(K);
end

if mod(K,2) == 0
    K = K + 1;
end

for t = 1:Ntest
    TestC{t} = test(t,:);
    Dist(t,:) = sqrt(sum((ones(Ntrain,1) * TestC{t} - train).^2,2)).';
end

[Y,I] = sort(Dist,2);

I = I(:,1:K);
Y = Y(:,1:K) + eps*rand(1);

for c = 1:Nclass
    for t = 1:Ntest
        for p = 1:K
            if find(I(t,p) == find(label == labelnames(c)))
                J(t,p) = labelnames(c);
            end
        end
    end
end

for c = 1:Nclass
    for t = 1:Ntest
        classscore(t,c) = length(find( J(t,:) == labelnames(c) ));
    end
end

Probability = classscore(:,wantid) ./ (sum(classscore,2));



