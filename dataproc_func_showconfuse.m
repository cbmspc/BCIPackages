function showconfuse (pconfuse)

me = mean(pconfuse,3) * 100;
sd = std(pconfuse,[],3) * 100;

Form = '% 5.1f';

for i = 1:size(me,1)
    for j = 1:size(me,2)
        fprintf('%.1f%% ± %.1f%% \t', me(i,j), sd(i,j));
    end
    fprintf('\n');
end

% for i = 1:size(me,1)
%     for j = 1:size(me,2)
%         cc{i,j} = [num2str(me(i,j),Form) '% ± ' num2str(sd(i,j),Form) '%'];
%     end
% end
% 
% cc
% 
