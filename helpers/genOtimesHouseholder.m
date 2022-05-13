%%  GENOTIMESHOUSEHOLDER Generates two cells comprising of a set of complex Hermitian or positive semidefinite matrices of the respective input dimensionalities
%   This function has five required input arguments:
%       FILENAME: file location where the cells are stored
%       DIM_A: integer that describes the dimension of the matrix for Alice
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       NUM_SUMMANDS: number of matrices for Alice and Bob each
%       OPT: "Pos" if positive semidefinite matrix is required
%            "Herm" if Hermitian matrix is required
%
%   genOtimesHouseholder(FILENAME, DIM_A, DIM_B, NUM_SUMMANDS, OPT) saves
%   NUM_SUMMANDS number of Hermitian or positive semidefinite matrices for
%   Alice and Bob respectively in FILENAME
%   Each matrix for Alice is of size DIM_A x DIM_A
%   Each matrix for Bob is of size DIM_B x DIM_B
%   FILENAME also contains the time required to generate each of the NUM_SUMMANDS matrices 
%   
%   URL: https://ankith-mohan.github.io/SEP/helpers/genOtimesHouseholder.html
%
%   requires: genHouseholder.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function genOtimesHouseholder(filename, dim_A, dim_B, num_summands, opt)
    K_list = cell(1, num_summands);
    L_list = cell(1, num_summands);
    
    K_time = zeros(1, num_summands);
    L_time = zeros(1, num_summands);
    
    h = waitbar(0, "Generating matrices");
    % Goal: \sum_{i=1}^{num_summands} K_i \otimes L_i
    for i = 1:num_summands
        
        waitbar(num_summands/i, h, sprintf("Summand %d: K", i));
        ki_start = tic;
        K_list{i} = genHouseholder(dim_A, opt);    
        K_time(i) = toc(ki_start);
        
        waitbar(num_summands/i, h, sprintf("Summand %d: L", i));
        li_start = tic;
        L_list{i} = genHouseholder(dim_B, opt);
        L_time(i) = toc(li_start);
    end
    
    save(filename, 'K_list', 'L_list', 'K_time', 'L_time', '-v7.3');
    close(h);
end