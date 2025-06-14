function [sphand, underhand, rsphand, handles_left_subplots_2dmatrix, handle_backgroundlayer, handles_right_subplots_2dmatrix] = subplotcompact(VN, HN, Opt)

HLS = 0.05;
HRS = 0.05;
HIS = 0.00;
VUS = 0.05;
VLS = 0.05;
VIS = 0.00;

if exist('Opt','var') && isstruct(Opt)
    if validate1(Opt,'LeftMargin')
        HLS = Opt.LeftMargin;
    end
    if validate1(Opt,'RightMargin')
        HRS = Opt.RightMargin;
    end
    if validate1(Opt,'HorizontalSeparation')
        HIS = Opt.HorizontalSeparation;
    end
    if validate1(Opt,'TopMargin')
        VUS = Opt.TopMargin;
    end
    if validate1(Opt,'BottomMargin')
        VLS = Opt.BottomMargin;
    end
    if validate1(Opt,'VerticalSeparation')
        VIS = Opt.VerticalSeparation;
    end
end

Hsize = (1.00 - HLS - HRS - (HN-1)*HIS) / HN;
Vsize = (1.00 - VLS - VUS - (VN-1)*VIS) / VN;

LeftAxesColor = 'none';
RightAxesColor = 'none';
UnderAxesColor = [1 1 1];

Posx = zeros(1,HN);
for i = 1:HN
    Posx(i) = HLS + HIS*(i-1) + Hsize*(i-1);
end
Posy = zeros(1,VN);
for i = 1:VN
    Posy(VN+1-i) = (VLS + VIS*(i-1) + Vsize*(i-1));
end

underhand = axes('Position', [HLS/2, VLS/2, 1-(HLS+HRS)/1.5, 1-(VLS+VUS)/1.5]);
set(underhand, 'XTick', [], 'YTick', [], 'Color', UnderAxesColor, 'XLim', [HLS/2, HLS/2+1-(HLS+HRS)/1.5], 'YLim', [VLS/2, VLS/2+1-(VLS+VUS)/1.5]);
handle_backgroundlayer = underhand;

if nargout >= 3
    rsphand = zeros(1,VN*HN);
    handles_right_subplots_2dmatrix = zeros(VN,HN);
    for k = reshape(flipud(reshape(1:VN*HN, HN, [])), 1, [])
        i = mod(k-1,HN)+1;
        j = floor((k-1)/HN)+1;
        rsphand(k) = axes('Position', [Posx(i), Posy(j), Hsize, Vsize], 'YAxisLocation', 'right', 'Color', RightAxesColor);
        handles_right_subplots_2dmatrix(j,i) = rsphand(k);
    end
end

sphand = zeros(1,VN*HN);
handles_left_subplots_2dmatrix = zeros(VN,HN);
for k = reshape(flipud(reshape(1:VN*HN, HN, [])), 1, [])
    i = mod(k-1,HN)+1;
    j = floor((k-1)/HN)+1;
    sphand(k) = axes('Position', [Posx(i), Posy(j), Hsize, Vsize], 'Color', LeftAxesColor);
    handles_left_subplots_2dmatrix(j,i) = sphand(k);
end



function v = validate1(Opt, fieldstring)
v = false;
if isfield(Opt,fieldstring) && isscalar(Opt.(fieldstring)) && isnumeric(Opt.(fieldstring)) && Opt.(fieldstring) >= 0
    v = true;
end

