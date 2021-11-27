
function [y1, y2]=Crossover(x1,x2)

%     alpha=rand(size(x1));
%     
%     y1=round(alpha.*x1+(1-alpha).*x2);
%     y2=round(alpha.*x2+(1-alpha).*x1);
    nvar = numel(x1);
    cut = randsample(1:nvar,1);
    y1 = [x1(1:cut),x2(cut+1:end)];
    y2 = [x2(1:cut),x1(cut+1:end)];
end