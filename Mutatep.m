% 
% 
 function y=Mutatep(x1,x2,x3,mu)
% 
%     nVar=numel(x);
%    
%     nMu=ceil(mu*nVar);
% 
%     j=randsample(nVar,nMu);        
    y=round(x1+mu*(x2-x3)/2);
   %y=intersect(union(x1,[mu]), setdiff(x2,x3));
 
   % y(j) = ~x(j);%+sigma.*randn(size(j));

end

% 
% function [ z ] = diffxy( x,y)
% % This function to find the difference in x not in y
% for i=1:length(x)
%     if x(i)~=y(i)
%         z(i)=x(i);
%     else
%         z(i)=0;
%     end
% end
% end
% y=or(x1,mu*(diffxy(x2,x3)));
%  end