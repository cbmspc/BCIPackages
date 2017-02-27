function ChanNumeric = chan2numeric (ChanNames)
if ~iscell(ChanNames)
    ChanNames = {ChanNames};
end
b = regexp(ChanNames, '^(?:(?:Grid)|(?:MG)|(?:G))(\d+)', 'tokens', 'once');
ChanNumeric = zeros(1,length(ChanNames));
for i = 1:length(b);
    if ~isempty(b{i})
        ChanNumeric(i) = str2double(b{i}{1});
    end
end

