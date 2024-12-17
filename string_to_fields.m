% Turns A=value, B=value, C=value, D=value, ... to a structure with fields A, B, C, D, ...
function fields = string_to_fields (str)
fields = struct;
b = regexp(str, '\w+=');
if ~isempty(b)
    b = [b length(str)+1];
    for i = 1:length(b)-1
        substr = str(b(i):b(i+1)-1);
        b2 = regexp(substr, '^([A-Za-z]\w+)=(.*)[,;]\s*$', 'tokens', 'once');
        if ~isempty(b2)
            fields.(b2{1}) = str2num(b2{2}); %#ok<ST2NM>
        end
    end
end

