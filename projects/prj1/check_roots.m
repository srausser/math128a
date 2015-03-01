% This function checks all the permutations of rts and exact, and outputs
% the smallest errors from all permutations.
%
% INPUT: rts = computed roots from your function
%        exact = 'exact' roots from the matlab ROOTS function
% OUTPUT: errors = smallest absolute error of any permutation of rts and
%                  exact
%         rts = same elements as rts, but reordered
%         exact = same elements as exact, but reordered
% USAGE EXAMPLE:
% rand_coeff = randn(1,4);
% rts = cubicx(rand_coeff);
% exact = roots(rand_coeff);
% err = check_roots(rts,exact);
function [errors, rts, exact, flg] = check_roots(rts, exact)

flg = 0;
if length(rts) ~= length(exact)
    flg = 1;
    errors = 0;
    return;
elseif sum(isnan(rts)) > 0
    flg = 2;
    errors = 0;
    return;
end

rts_perms = perms(rts);
exact_perms = perms(exact);

m = size(rts_perms,1);
n = size(exact_perms,1);

min_err = Inf;

for i=1:m
    for j=1:n
        vec1 = rts_perms(i,:);
        vec2 = exact_perms(j,:);
        if norm(vec1-vec2) < min_err
            min_err = norm(vec1-vec2);
            errors = abs(vec1-vec2)./max(vec2);
            rts = vec1;
            exact = vec2;
        end
    end
end

end