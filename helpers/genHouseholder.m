function M = genHouseholder(dim, opt)
    % GENPOSHOUSEHOLDER Generate a complex Hermitian matrix of desired dimensionality 
    %   H = GENPOSHOUSEHOLDER(dim) generates a dim x dim complex Hermitian matrix
    
    if strcmp(opt, "Pos")
        % Eigenvalues are ``dim" numbers sampled uniformly from [0,1] for PSD
        lambda = rand(1, dim);
    elseif strcmp(opt, "Herm")
        a = -1;
        b = 1;
        % Eigenvalues are ``dim" numbers sampled uniformly from [-1,1] for Herm
        lambda = (b - a) .* rand(1, dim) + a;
    else
        error("Expected 'Herm' or 'Pos' for `opt`, got %s instead.", opt);
    end
    Lambda = diag(lambda); % Matrix of eigenvalues
    %%%%%%%%%%%%%%%%%%%% BEWARE OF THIS EXPRESSION %%%%%%%%%%%%%%%%%%%%%%%%
    M = complex(rand(dim), rand(dim)); % Create a complex matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Q, ~] = qr(M); % Use QR decomposition to make the columns of M orthonormal
    M = Q * Lambda * Q'; % This should make M have a good spread of eigenvalues
    M = 0.5 * (M + M'); % This should guarantee Hermiticity (especially on the principal diagonal)
    if strcmp(opt, "Pos")
        assert(IsPSD(M) == 1, "Resultant matrix is not PSD"); % Are you really sure its PSD?
    else
        assert(ishermitian(M), "Resultant matrix is not Herm"); % Are you really sure its Herm?
    end
end