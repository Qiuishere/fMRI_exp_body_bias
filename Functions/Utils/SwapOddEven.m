function Swapped = SwapOddEven(In,Dim)

% swap odd and even rows/columns
% i.e. 1 goes to 2 and 2 to 1, 3 to 4 and 4 to 3, etc.
% In: array to rearrange
% Dim: dimension along which to rearrange (1 or 2)

N = size(In,Dim);
K = rem(N,2);

if K==1
    warning('length is odd! Last row/column will be left unaltered');
end

ind = 1;
Ordr = nan(1,N);

for i = 1:ceil(N/2)
    Ordr(ind) = ind + 1;
    Ordr(ind + 1) = ind;
    ind = ind + 2;
end

switch Dim
    case 1
        Swapped = In(Ordr,:);
    case 2
        Swapped = In(:,Ordr);
end

end
