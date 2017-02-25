function [trfmat, Ftrain, Flabel, SuffStat, latent] = cpca_ida(TrainData, TrainLabels, Fdim)
cmat = dataproc_func_cpca(TrainData,TrainLabels,[]);
for i = 1:length(cmat)
    amat{i} = dataproc_func_ida(TrainData*cmat{i},TrainLabels,Fdim);
    trfmat{i} = cmat{i}*amat{i};
    Ftrain{i} = TrainData*trfmat{i};
end
Flabel = TrainLabels;

if nargout >= 4
    SuffStat.classes = unique(Flabel);
    for s = 1:length(cmat)
        for c = 1:length(SuffStat.classes)
            SuffStat.stats{s,c,1} = mean(Ftrain{s}(find(Flabel==SuffStat.classes(c)),:),1);
            SuffStat.stats{s,c,2} = cov(Ftrain{s}(find(Flabel==SuffStat.classes(c)),:));
        end
    end
end

latent = [];
