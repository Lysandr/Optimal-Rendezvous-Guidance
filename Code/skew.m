function skewed = skew(mat)

if numel(mat) >  3
    error('Vector size can only be 3');
end

% provides the skew symmetric form of a matrix
%
    skewed = [0 -mat(3) mat(2); mat(3) 0 -mat(1); -mat(2) mat(1) 0];


end