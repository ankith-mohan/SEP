%%  EXECUTEEXPT Creates the necessary folder structure for SETUPEXPT
%   This function has sixteen required input arguments:
%       DIM_A: integer that describes the dimension of the matrix for Alice
%       DIM_B: integer that describes the dimension of the matrix for Bob
%       K_SE: integer that describes the maximum number of symmetric extensions
%       K_BSE: integer that describes the maximum number of bosonic
%       extensions
%       N_SEESAW: maximum number of see-saw steps
%       N_RAND: number of random starting points
%       SUM: 1 if GENOTIMESHOUSEHOLDER is used to generate the cells of
%       matrices
%            0 otherwise
%       NUM_SUMMAND: number of matrices that have to be summed over
%       N_OPS: number of instances of the input matrix to compute the
%       bounds for
%       THRESHOLD: threshold to determine convergence
%       OPT: "sigma_A" if starting matrix is for Alice's subsystem
%            "sigma_B" is starting matrix if for Bob's subsystem
%       PROJ: 0 to use GENHOUSEHOLDER
%             1 to use GENHERM if OPT == "Herm"
%                      GENPOS if OPT == "Pos"
%       DONT_DO: list that contains what functions must NOT be called
%                Upper bounds:
%                "PPT": computebeta_PPT if SUM == 0 else compute_PPT_sum
%                "realignment": computebeta_r if SUM == 0 else
%                computebeta_r_sum
%
%                Symmetric Extensions:
%                "SymExtOnly": computebeta_k(PPT=0) if SUM == 0 else
%                computebeta_k_sum(PPT=0)
%                "SymExtOnlyrlg": computebeta_kr(PPT=0) if SUM == 0 else
%                computebeta_kr_sum(PPT=0)
%                "SymExt": computebeta_k(PPT=1) if SUM == 0 else
%                computebeta_k_sum(PPT=1)
%                "SymExtrlg": computebeta_kr(PPT=1) if SUM == 0 else
%                computebeta_kr_sum(PPT=1)
%               
%                Bosonic Extensions:
%                "BosSymExtOnly": computebeta_prime_k(PPT=0) if SUM == 0
%                else computebeta_prime_k_sum(PPT=0)
%                "BosSymExtOnlyrlg": computebeta_prime_kr(PPT=0) if SUM ==
%                0 else computebeta_prime_kr_sum(PPT=0)
%                "BosSymExt": computebeta_prime_k(PPT=1) if SUM == 0 else
%                computebeta_prime_k_sum(PPT=1)
%                "BosSymExtrlg": computebeta_prime_kr(PPT=1) if SUM == 0
%                else computebeta_prime_kr_sum(PPT=1)
%
%                Lower bounds:
%                Maximally mixed state:
%                "MM_simple": computegamma_part(SIGMA=eye(dim_A)/dim_A) if
%                SUM == 0 else computegamma_part_sum(SIGMA=eye(dim_A)/dim_A)
%                "MM_rev_1": computegamma_part_rev_1(SIGMA=eye(dim_A)/dim_A) if
%                SUM == 0 else computegamma_part_rev_1_sum(SIGMA=eye(dim_A)/dim_A)
%                "MM_rev_2": computegamma_part_rev_2(SIGMA=eye(dim_A)/dim_A) if
%                SUM == 0 else computegamma_part_rev_2_sum(SIGMA=eye(dim_A)/dim_A)
%
%                State of uniform superposition:
%                "US_simple": computegamma_part(SIGMA=ones(dim_A)/dim_A) if
%                SUM == 0 else computegamma_part_sum(SIGMA=ones(dim_A)/dim_A)
%                "US_rev_1": computegamma_part_rev_1(SIGMA=ones(dim_A)/dim_A) if
%                SUM == 0 else computegamma_part_rev_1_sum(SIGMA=ones(dim_A)/dim_A)
%                "US_rev_2": computegamma_part_rev_2(SIGMA=ones(dim_A)/dim_A) if
%                SUM == 0 else computegamma_part_rev_2_sum(SIGMA=ones(dim_A)/dim_A)
%
%                Random starting points:
%                "rand_simple": computegamma_rand(SS_TYPE="simple") if SUM 
%                == 0 else computegamma_rand_sum(SS_TYPE="simple")
%                "rand_rev_1": computegamma_rand(SS_TYPE="rev_1") if SUM ==
%                0 else computegamma_rand_sum(SS_TYPE="rev_1")
%                "rand_rev_2": computegamma_rand(SS_TYPE="rev_2") if SUM ==
%                0 else computegamma_rand_sum(SS_TYPE="rev_2")
%
%       SHOULD_SAVE: 1 if the density matrices returned by each of the
%       functions need to be stored in respective file in FOLDERNAME
%       FOLDERNAME: folder where the file with the results is to be created
%       PLUS_TWO_RAND: 2 if both the maximally mixed state and the state of
%       uniform superposition must be added to the list of starting points
%                      0 otherwise
%
%   executeExpt(DIM_A, DIM_B, K_SE, K_BSE, N_SEESAW, N_RAND, SUM, NUM_SUMMAND, 
%   N_OPS, THRESHOLD, OPT, PROJ, DONT_DO, SHOULD_SAVE, FOLDERNAME, PLUS_TWO_RAND)
%   writes the computed values in 
%   ./FOLDERNAME/Alice=DIM_A,Bob=DIM_B,OPT.mat if SUM == 0 else
%   ./FOLDERNAME/Alice=DIM_A,Bob=DIM_B/op_j/Alice=DIM_A,Bob=DIM_B,OPT.mat
%   for j = 1:N_OPS
%   
%   URL: https://ankith-mohan.github.io/SEP/main/executeExpt.html
%
%   requires: setupExpt.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022  


function executeExpt(dim_A, dim_B, K_se, K_bse, N_seeSaw, N_rand, sum, ...
                     num_summand, N_ops, threshold, opt, proj, dont_do, ...
                     should_save, foldername, plus_two_rand)
    name = sprintf("Alice=%d,Bob=%d", dim_A, dim_B);
    suffix = sprintf("%s,%s", name, opt);
    if sum == 1
        reqd_folder = sprintf("./%s/%s", foldername, name);
        if not(isfolder(reqd_folder))
            mkdir(reqd_folder);
        end
        
        for j = 1:N_ops
            reqd_subfolder = sprintf("%s/op_%d", reqd_folder, j);
            if not(isfolder(reqd_subfolder))
                mkdir(reqd_subfolder);
            end
            
            filename = sprintf("%s/%s", reqd_subfolder, suffix);
            setupExpt(dim_A, dim_B, K_se, K_bse, N_seeSaw, N_rand, sum, ...
                      num_summand, N_ops, threshold, opt, filename, ...
                      proj, dont_do, should_save, plus_two_rand);               
        end
    else
        filename = sprintf("./%s/%s", foldername, suffix);
        setupExpt(dim_A, dim_B, K_se, K_bse, N_seeSaw, N_rand, sum, ...
              num_summand, N_ops, threshold, opt, filename, proj, ...
              dont_do, should_save, plus_two_rand);
    end
end