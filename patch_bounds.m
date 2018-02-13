% This function creates patches in the specified Bounds
% Use this function before plotting time series EEG

% You can color each bound in different color, just put them in separate
% rows.

function phand = patch_bounds (AxesHand, Bounds, Ylim, FaceColor, EdgeColor, FaceAlpha, EdgeAlpha)
if isempty(AxesHand)
    AxesHand = gca;
end

if ~exist('Ylim','var') || isempty(Ylim) || size(Ylim,2) ~= 2
    Ylim = [0 1];
end

if ~exist('FaceColor','var') || isempty(FaceColor)
    FaceColor = [0.9 0.9 1.0];
end

if ~exist('EdgeColor','var') || isempty(EdgeColor)
    EdgeColor = FaceColor;
end

if ~exist('FaceAlpha','var') || isempty(FaceAlpha)
    FaceAlpha = 1;
end

if ~exist('EdgeAlpha','var') || isempty(EdgeAlpha)
    EdgeAlpha = 1;
end

hold(AxesHand,'on');

if size(FaceColor,1) ~= size(Bounds,1)
    FaceColor = ones(size(Bounds,1),1) * FaceColor(1,:);
end

if size(EdgeColor,1) ~= size(Bounds,1)
    EdgeColor = ones(size(Bounds,1),1) * EdgeColor(1,:);
end

if size(FaceAlpha,1) ~= size(Bounds,1)
    FaceAlpha = ones(size(Bounds,1),1) * FaceAlpha(1,:);
end

if size(EdgeAlpha,1) ~= size(Bounds,1)
    EdgeAlpha = ones(size(Bounds,1),1) * EdgeAlpha(1,:);
end

if size(Ylim,1) ~= size(Bounds,1)
    for i = 2:size(Bounds,1)
        Ylim(i,:) = Ylim(1,:);
    end
end

phand = zeros(1,size(Bounds,1));
axes(AxesHand);
for i = 1:size(Bounds,1)
    phand(i) = patch([Bounds(i,1) Bounds(i,2) Bounds(i,2) Bounds(i,1)], [Ylim(i,1) Ylim(i,1) Ylim(i,2) Ylim(i,2)], FaceColor(i,:), 'EdgeColor', EdgeColor(i,:), 'FaceAlpha', FaceAlpha(i,1), 'EdgeAlpha', EdgeAlpha(i,1));
end

