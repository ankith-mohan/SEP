%%  GENHERM  Generates a complex Hermitian matrix of the input dimensionality
%   This function has one required input argument:
%       DIM: integer that describes the local dimension of the matrix
%
%   H = genHerm(DIM) is a Hermitian matrix of size dim x dim
%


function H = genHerm(dim)
    % GENHERM Generate a complex Hermitian matrix of desired dimensionality 
    %   H = GENHERM(dim) generates a dim x dim complex Hermitian matrix
    
    M = complex(rand(dim), rand(dim)); % Create a complex matrix
    H = 0.5 * (M + M'); % Make it Hermitian
    assert(ishermitian(H), "H is not Hermitian"); % Are you really sure its Hermitian?
end