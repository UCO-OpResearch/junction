function [Pred Rank] = InitializeRootedTree(n)
% This function creates a predecessor list and a rank list for "n" nodes in
% a rooted tree forest. At this point all nodes are in their own tree so
% pred(i)=i and Rank(i)=0

for i = 1:n
    Pred(i) = i;
    Rank(i) = 0;
end
    