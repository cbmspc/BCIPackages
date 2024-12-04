% Opposite of getbounds
%
% Example usages with simple bounds:
%
%   Create a vector of length 30, where indices 1 thru 10 and 21 thru 30
%   have the value True, and the rest have value False:
%     putbounds([1 10; 21 30])
%     The first input argument is required.
%
%   If the first input argument has three columns instead of two, instead
%   of only putting True and False in the output vector, put the values
%   specified in the third column (41 and 52 in this example):
%     putbounds([1 10 41; 21 30 52])
%
%   2nd input [optional]: With length 40:
%     putbounds([1 10; 21 30], 40)
%     If left empty, the default is the minimum length necessary to hold
%     both the bounds and prepads (prepads are defined by the 3rd input)
%
%   3rd input [optional]: Still keeping length to 40, insert 7 default
%   values to the front, i.e. prepadding.
%     putbounds([1 10; 21 30], 40, 7)
%     If left empty, do not insert anything in front
%
%   4th input [optional]: Instead of defaulting to False, the default value
%   is now 42 (this automatically upgrades the output vector to the data
%   type necessary to hold this kind of value):
%     putbounds([1 10; 21 30], 40, 7, 42)
%     If left empty, the default default value is False or 0 for prepadding
%     and in between specified bounds
%


function Signal = putbounds (Bounds, SigLen, PrePadThisManySamples, DefaultValue)
if ~exist('PrePadThisManySamples','var') || isempty(PrePadThisManySamples) || PrePadThisManySamples < 0
    PrePadThisManySamples = 0;
end
if ~exist('DefaultValue','var') || isempty(DefaultValue) || ~isscalar(DefaultValue)
    DefaultValue = false;
end
if ~exist('SigLen','var') || isempty(SigLen) || ~isscalar(SigLen) || ~isnumeric(SigLen)
    SigLen = max(Bounds(:,2)) + PrePadThisManySamples;
else
    SigLen = round(SigLen);
end

if SigLen <= 0
    Signal = repmat(DefaultValue, 0);
    return
end

isTri = false;
if size(Bounds,2) >= 3
    if DefaultValue ~= 0
        Signal = repmat(DefaultValue,[SigLen 1]);
    else
        Signal = zeros(SigLen,1);
    end
    isTri = true;
else
    Signal = false(SigLen,1);
end


for i = 1:size(Bounds,1)
    a = round(Bounds(i,1));
    b = round(Bounds(i,2));
    c = true;
    if isTri
        c = Bounds(i,3);
    end
    if a < 1
        a = 1;
    end
    if a > SigLen
        a = SigLen;
    end
    Signal(a:b) = c;
end

if PrePadThisManySamples > 0
    if exist('PrePadWithThisValue','var') && isscalar(DefaultValue) && (isnumeric(DefaultValue) || islogical(DefaultValue))
        c = DefaultValue;
    else
        c = false;
    end
    Signal = [repmat(c, [round(PrePadThisManySamples) 1]); Signal];
end

if length(Signal) > SigLen
    Signal = Signal(1:SigLen);
end
