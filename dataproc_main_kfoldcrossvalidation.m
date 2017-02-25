% Data = Time domain data (Time, Chan)
% Labels = Discrete or continuous labeling of Data
% K = K fold
% Func = Function handle that does the actual validation
function Results = dataproc_main_kfoldcrossvalidation (Data, Labels, K, Func, Opts)

Ndatatime = size(Data,1);
Nlabelstime = size(Labels,1);

DataTimeDivs = (0:Ndatatime/K:Ndatatime)';
DataBounds = round([DataTimeDivs(1:end-1)+1 DataTimeDivs(2:end)]);

LabelsTimeDivs = (0:Nlabelstime/K:Nlabelstime)';
LabelsBounds = round([LabelsTimeDivs(1:end-1)+1 LabelsTimeDivs(2:end)]);

Results = cell(K,1);

for k = K:-1:1
    % k being the test set
    DataTimePointsTest = DataBounds(k,1):DataBounds(k,2);
    DataTimePointsTrain = setdiff(1:Ndatatime, DataTimePointsTest);
    
    LabelsTimePointsTest = LabelsBounds(k,1):LabelsBounds(k,2);
    LabelsTimePointsTrain = setdiff(1:Nlabelstime, LabelsTimePointsTest);
    
    Results{k} = Func(Data(DataTimePointsTrain,:), Labels(LabelsTimePointsTrain,:), Data(DataTimePointsTest,:), Labels(LabelsTimePointsTest,:), Opts);
    
    
end
