function H = genHerm(dim)
    % GENHERM Generate a complex Hermitian matrix of desired dimensionality 
    %   H = GENHERM(dim) generates a dim x dim complex Hermitian matrix
    
    H = complex(rand(dim), rand(dim)); % Create a complex matrix
    H = 0.5 * (H + H'); % Make it Hermitian
    assert(ishermitian(H) == 1, "H is not Hermitian"); % Are you really sure its Hermitian?
end