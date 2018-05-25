function y = bridge_nans (x, method, extrap)
if ~exist('method','var') || isempty(method)
    method = 'linear';
end
if ~exist('extrap','var')
    extrap = 'extrap';
end

y = x;

if ~isempty(extrap)
    for ch = 1:size(x,2)
        y(:,ch) = interp1(find(~isnan(x(:,ch))), x(~isnan(x(:,ch)),ch), (1:size(x,1)).', method, extrap);
        
        %     nn = isnan(x(:,ch));
        %     nnbb = getdigitalbounds(nn);
        %     if ~isempty(nnbb)
        %         for i = 1:size(nnbb,1)
        %             mark = nnbb(i,:);
        %             XX = [max(1,mark(1)-3):mark(1)-1, mark(2)+1:min(size(x,1),mark(2)+3)];
        %             y(mark(1):mark(2),ch) = interp1(XX, x(XX, ch), [mark(1):mark(2)].', 'pchip');
        %         end
        %     end
    end
else
    for ch = 1:size(x,2)
        k = find(~isnan(x(:,ch)));
        if ~isempty(k)
            y(:,ch) = interp1(k, x(~isnan(x(:,ch)),ch), (1:size(x,1)).', method);
            y(1:k(1)-1,ch) = NaN;
            y(k(end)+1:end,ch) = NaN;
        end
    end
end
return
