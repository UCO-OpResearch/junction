function [Pred Root] = FindRoot(Pred, Node)
% This function accepts A Predecessor list for a rooted forest and a Node.
% The function then finds the root of the tree containing "Node" and, in
% the process, uses path halving to shorten subsequent searches.

Root = Node; % Set the root equal to the start Node
while Pred(Root) ~= Root % as long as the predecessor of the root is not the root itself
    OldNode = Root; % Store the Node we're currently at
    Root = Pred(Root); % Move one step up the tree
    Pred(OldNode) = Pred(Root); % Set the predecessor of the stored node to the node two levels above it
    Root = Pred(Root); % Move another step up the tree
end
Root = Pred(Root); % Once the "while" loop exits we could still be one level below the root so we'll move up one more step.
% NOTE: we may make up to two more steps than we need to but this won't
% be a problem since Pred(Root) = Root so we'll just cycle for a few steps.
% Also the additional steps caused by this are well worth the savings from
% the path halving.