function Ord = SortToMatch(X,Y)
% returns ordering needed to match vector X to vector Y.

assert(length(X)==length(Y), 'Length of two vectors must be same!');

Ord = zeros(size(X));

for i = 1:numel(X)
   Ord(i) = find(X==Y(i),1,'first'); % find where element Y(i) is in X
   X(find(X==Y(i),1,'first')) = NaN; % remove element from X
end


end