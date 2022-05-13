%%  GENHOUSEHOLDER  Generates a complex Hermitian or positive semidefinite matrix of the input dimensionality
%   This function has two required arguments:
%       DIM: integer that describes the local dimension of the matrix
%       OPT: "Pos" if positive semidefinite matrix is required
%            "Herm" if Hermitian matrix is required
%
%   M = genHouseholder(DIM, OPT) is a Hermitian or positive semidefinite
%   matrix of size DIM x DIM
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/genHouseholder.html
%
%   requires: IsPSD (http://qetlab.com/IsPSD)
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function M = genHouseholder(dim, opt)
    if strcmp(opt, "Pos")
        % Eigenvalues are DIM numbers sampled uniformly from [0,1] for PSD
        lambda = rand(1, dim);
    elseif strcmp(opt, "Herm")
        a = -1;
        b = 1;
        % Eigenvalues are DIM numbers sampled uniformly from [-1,1] for Herm
        lambda = (b - a) .* rand(1, dim) + a;
    else
        error("Expected 'Herm' or 'Pos' for `opt`, got %s instead.", opt);
    end
    % Diagonal matrix of eigenvalues
    Lambda = diag(lambda);
    % Create a complex matrix
    M = complex(rand(dim), rand(dim));
    % Use QR decomposition to make the columns of M orthonormal
    [Q, ~] = qr(M);
    % This should make M have a good spread of eigenvalues
    M = Q * Lambda * Q';
    % This should guarantee Hermiticity (especially on the principal diagonal)
    M = 0.5 * (M + M');
    if strcmp(opt, "Pos")
        % This should make it PSD
        M = M' * M;
        % Are you really sure its PSD?
        assert(IsPSD(M) == 1, "Resultant matrix is not PSD");
    else
        % Are you really sure its Herm?
        assert(ishermitian(M), "Resultant matrix is not Herm");
    end
end