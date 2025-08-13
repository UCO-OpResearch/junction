function [Pred Rank] = UnionbyRank(Pred, Rank, Tree1, Tree2)
% This function accepts a Predecessor and Rank function describing a rooted
% tree forest and two trees to be merged. It then merges the trees based on
% their ranks and updates the Pred and Rank lists

if Rank(Tree1) > Rank(Tree2)
    Pred(Tree2) = Tree1; % Merge the smaller tree into the larger tree
elseif Rank(Tree1) < Rank(Tree2)
    Pred(Tree1) = Tree2; % Merge the smaller tree into the larger tree
else
    Pred(Tree1) = Tree2; % Merge the trees arbitrarily
    Rank(Tree2) = Rank(Tree2) + 1; % Increase the rank of the new tree by 1
end