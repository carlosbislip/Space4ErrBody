%--------------------------------------------------------------------------
% This function identifies the pareto frontier of a set of points (this
% function consider smaller values are more desirable)
%--------------------------------------------------------------------------
% Input: input, a matrix, each row correspondes to a point, each column
% correspond to a dimension
%--------------------------------------------------------------------------
% Outputs:
% (1) membership: a logical array, have same number of rows as input
% matrix, 1 indicate the corresponding point in input matrix is a member of
% pareto frontier, 0 otherwise
% (2) member_value: matrix, contain point(s) on the pareto frontier.
%--------------------------------------------------------------------------
% Example:
%   x=rand(100,2);
%   [membership,member_value]=find_pareto_frontier(x);
%   scatter(x(:,1),x(:,2));
%   hold on;
%   scatter(member_value(:,1),member_value(:,2),'r');
%   legend({'Data','Pareto Frontier'})
%
%--------------------------------------------------------------------------
function [membership, member_value]=find_pareto_frontier(input)
out=[];

data=unique(input,'rows');
for i = 1:size(data,1)
    
    c_data = repmat(data(i,:),size(data,1),1);
    t_data = data;
    t_data(i,:) = Inf(1,size(data,2));
    smaller_idx = c_data>=t_data;
    
    idx=sum(smaller_idx,2)==size(data,2);
    if ~nnz(idx)
        out(end+1,:)=data(i,:);
    end
end
membership = ismember(input,out,'rows');
member_value = out;
end