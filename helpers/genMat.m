%% GENMAT Generates a complex Hermitian or positive semidefinite matrix of the input dimensionality normalized by its trace norm
%   This function has three required input arguments:
%       DIM: integer that describes the local dimension of the matrix
%       OPT: "Pos" if positive semidefinite matrix is required
%            "Herm" if Hermitian matrix is required
%       PROJ: 0 to use GENHOUSEHOLDER
%             1 to use GENHERM if OPT == "Herm" or GENPOS if OPT == "Pos"
%
%   M = genMat(DIM, OPT, PROJ) is a Hermitian or positive semidefinite
%   matrix of size DIM x DIM
%
%   URL: https://ankith-mohan.github.io/SEP/helpers/genMat.html
%
%   requires: genHerm.m, genPos.m, genHouseholder.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function M = genMat(dim, opt, proj)
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
        % Trace norm
        M_op = max(abs(eig(M)));
        % Normalize by trace norm
        M = M/M_op;
    end
end