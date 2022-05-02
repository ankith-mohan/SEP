function M = genMat(dim, opt, proj)
    % GENMAT Generate either a dim x dim complex Hermitian or a dim x dim
    % complex positive semidefinite matrix normalized by its trace norm
    % M = GENMAT(dim, opt, proj) generates a dim x dim complex Hermitian matrix
    % if `opt` is "Herm", or a dim x dim complex positive semidefinite
    % matrix if `opt` is "Pos"
    % if `proj` is 1 then generate a Projection matrix
    
    if strcmp(opt, "Herm")
        if proj == 1
            M = genHerm(dim);
        else
            M = genHouseholder(dim, "Herm");
        end
    elseif strcmp(opt, "Pos")
        if proj == 1
            M = genPos(dim);
        else
            M = genHouseholder(dim, "Pos");
        end
    else
        error("Expected 'Herm' or 'Pos' for `opt`, got %s instead.", opt);
    end
    if proj == 1
        M_op = max(abs(eig(M))); % Trace norm
        M = M/M_op; % Normalize by trace norm
    end
end