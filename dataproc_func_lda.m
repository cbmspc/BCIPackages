
% Feature extraction based on Linear Discriminant Analysis
% TrainData is (obs x dim) (Each row is an observation)

% Feature dimension is limited to Nclass-1

function Fmat = featureextraction(TrainData, TrainLabels, MfeatureSpace)
%

if any(isnan(TrainData(:)))
    error('TrainData contains NaN');
end

if any(isnan(TrainLabels(:)))
    error('TrainLabels contains NaN');
end


if isempty(MfeatureSpace)
    warning('MfeatureSpace is blank (not set). Results may be unpredictable.');
end

BinData = double(TrainData).';
classes = unique(TrainLabels);
Nclass = length(classes);
Nparm = size(TrainData,2);

if MfeatureSpace > Nclass-1
    MfeatureSpace = Nclass-1;
end

for i = 1:Nclass
    Memship(find(TrainLabels == classes(i))) = i-1;
end

n = size(BinData,1);
m = MfeatureSpace;

M = zeros(n,1);
Sw = zeros(n,n);
Sb = zeros(n,n);
Si = zeros(Nparm,Nparm,Nclass);
for j = 1:Nclass
    ind = find(Memship == j-1);
    p_i(j) = length(ind)/length(Memship);
    Mi(:,j) = mean(BinData(:,ind),2);
    if length(ind) > 1
        Si(:,:,j) = cov(BinData(:,ind)',1); %20160205: Changed to pop cov
    else
        Si(:,:,j) = zeros(Nparm,Nparm);
    end
    %overall mean
    M = M + p_i(j) * Mi(:,j);

    %within class matrix
    Sw = Sw + p_i(j) * Si(:,:,j);

    %between class matrix
    Sb = Sb + p_i(j) * [Mi(:,j) * Mi(:,j)'];
end

Sb = Sb - M * M';
S = Sw + Sb;

%display options for eigs
OPTS.disp = 0;

try
    [V, D]= eigs(Sb,Sw,m,'LM',OPTS);
catch
    [V, D]= eig(Sb,Sw);
end
Fmat = V/norm(V);
