function [trfmat, Ftrain, Flabel, SuffStat, latent] = cpca_lda(TrainData, TrainLabels, Fdim)
cmat = dataproc_func_cpca(TrainData,TrainLabels,[]);
amat = cell(1,length(cmat));
trfmat = amat;
Ftrain = amat;
for i = 1:length(cmat)
    amat{i} = dataproc_func_lda(TrainData*cmat{i},TrainLabels,Fdim);
    trfmat{i} = cmat{i}*amat{i};
    Ftrain{i} = TrainData*trfmat{i};
end
Flabel = TrainLabels;

if nargout >= 4
    SuffStat.classes = unique(Flabel);
    for s = 1:length(cmat)
        for c = 1:length(SuffStat.classes)
            SuffStat.stats{s,c,1} = mean(Ftrain{s}(Flabel==SuffStat.classes(c),:),1);
            SuffStat.stats{s,c,2} = cov(Ftrain{s}(Flabel==SuffStat.classes(c),:));
        end
    end
end

% As CPCA is piecewise, latent is not applicable.
latent = [];