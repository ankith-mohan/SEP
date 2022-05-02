function P = genPos(dim)
    % GENPOS Generate a complex positive semidefinite matrix of desired dimensionality 
    %   P = GENPOS(dim) generates a dim x dim complex positive semidefinite matrix
    
    H = genHerm(dim); % Create a Hermitian matrix
    P = H' * H; % Make it PSD
    assert(IsPSD(P), "P is not PSD"); % Sure its PSD?
end