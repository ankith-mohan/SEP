function A = genMat(dim, opt)
    % GENMAT Generate either a dim x dim complex Hermitian or a dim x dim
    % complex positive semidefinite matrix normalized by its trace norm
    % M = GENMAT(dim, opt) generates a dim x dim complex Hermitian matrix
    % if `opt` is "Herm", or a dim x dim complex positive semidefinite
    % matrix if `opt` is "Pos"
    
    if strcmp(opt, "Herm")
        A = genHerm(dim);
    elseif strcmp(opt, "Pos")
        A = genPos(dim);
    else
        error("Expected 'Herm' or 'Pos' for `opt`, got %s instead.", opt);
    end
    A_op = max(abs(eig(A))); % Trace norm
    A = A/A_op; % Normalize by trace norm
end