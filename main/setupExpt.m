%%  SETUPEXPT   Computes the upper and lower bounds for N_OPS instances
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
%       NUM_SUMMANDS: number of matrices that have to be summed over
%       N_OPS: number of instances of the input matrix to compute the
%       bounds for
%       THRESHOLD: threshold to determine convergence
%       OPT: "sigma_A" if starting matrix is for Alice's subsystem
%            "sigma_B" is starting matrix if for Bob's subsystem
%       FILENAME: ".mat" file where all the results have to written
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
%       functions need to be stored in FILENAME
%                    0 otherwise
%       PLUS_TWO_RAND: 2 if both the maximally mixed state and the state of
%       uniform superposition must be added to the list of starting points
%                      0 otherwise
%
%   setupExpt(DIM_A, DIM_B, K_SE, K_BSE, N_SEESAW, N_RAND, SUM, NUM_SUMMANDS, 
%   N_OPS, THRESHOLD, OPT, FILENAME, PROJ, DONT_DO, SHOULD_SAVE, PLUS_TWO_RAND)
%   writes the computed values in FILENAME.mat
%   
%   URL: https://ankith-mohan.github.io/SEP/main/setupExpt.html
%
%   requires: genOtimesHouseholder.m, genMat.m, computebeta_PPT.m,
%   computebeta_PPT_sum.m, computebeta_r.m, computebeta_r_sum.m,
%   computebeta_k.m, computebeta_k_sum.m, computebeta_kr.m,
%   computebeta_kr_sum.m, computebeta_prime_k.m, computebeta_prime_k_sum.m,
%   computebeta_prime_kr.m, computebeta_prime_kr_sum.m,
%   computegamma_part.m, computegamma_part_sum.m,
%   computegamma_part_rev_1.m, computegamma_part_rev_1_sum.m,
%   computegamma_part_rev_2.m, computegamma_part_rev_2_sum.m,
%   computegamma_rand.m, computegamma_rand_sum.m
%   author: Ankith Mohan (ankithmo@vt.edu)
%   last updated: May 2, 2022


function setupExpt(dim_A, dim_B, K_se, K_bse, N_seeSaw, N_rand, ...
                   sum, num_summands, N_ops, threshold, opt, filename, ...
                   proj, dont_do, should_save, plus_two_rand)
    % SETUPEXPT
    
    if should_save == 1
        warning("`should_save` has been set to `true`. This operation is resource expensive!");
    end
    
    if sum == 1
        N_ops = 1;
    end
    
    % Get the list of all the variables in the filename (if it exists)
    var_list = {};
    filename_ext = sprintf("%s.mat", filename);
    if isfile(filename_ext)
        var_list = who('-file', filename_ext);
    end
    
    %% Generate operator
    if sum == 1
        if any(ismember('K_list', var_list))
            load(filename_ext, 'K_list', 'L_list');
        end
    else
        % Load the operators (if they exist)
        if any(ismember('Pis', var_list)) 
            load(filename_ext, 'Pis');
        else
            Pis = cell(1, N_ops);
        end
    end
    
    %% PPT
    % Load the PPT values (if they exist)
    % Else initialize only if computation required
    if any(ismember('beta_PPT', var_list))
        load(filename, 'beta_PPT', 'beta_PPT_time');
    elseif ~any(ismember('PPT', dont_do))
        if should_save == 1
            rho_PPT = cell(1, N_ops);
        end
        beta_PPT = zeros(1, N_ops);
        beta_PPT_time = zeros(1, N_ops);
    end
    
    %% Realignment
    if any(ismember('beta_r', var_list))
        load(filename, 'beta_r', 'beta_r_time');
    elseif ~any(ismember('realignment', dont_do))
        if should_save == 1
            rho_r = cell(1, N_ops);
        end
        beta_r = zeros(1, N_ops);
        beta_r_time = zeros(1, N_ops);
    end
    
    %% Symmetric Extensions
    %% Without PPT
    % Without realignment
    if any(ismember('beta_konly', var_list))
        load(filename, 'beta_konly', 'beta_konly_time');
    elseif ~any(ismember('SymExtOnly', dont_do))
        if should_save == 1
            rho_konly = cell(K_se, N_ops);
        end
        beta_konly = zeros(K_se, N_ops);
        beta_konly_time = zeros(K_se, N_ops);
    end

    % With realignment
    if any(ismember('beta_konlyr', var_list))
        load(filename, 'beta_konlyr', 'beta_konlyr_time');
    elseif ~any(ismember('SymExtOnlyrlg', dont_do))
        if should_save == 1
            rho_konlyr = cell(K_se, N_ops);
        end
        beta_konlyr = zeros(K_se, N_ops);
        beta_konlyr_time = zeros(K_se, N_ops);
    end

    %% With PPT
    % Without realignment
    if any(ismember('beta_k', var_list))
        load(filename, 'beta_k', 'beta_k_time');
    elseif ~any(ismember('SymExt', dont_do))
        if should_save == 1
            rho_k = cell(K_se, N_ops);
        end
        beta_k = zeros(K_se, N_ops);
        beta_k_time = zeros(K_se, N_ops);
    end

    % With realignment
    if any(ismember('beta_kr', var_list))
        load(filename, 'beta_kr', 'beta_kr_time');
    elseif ~any(ismember('SymExtrlg', dont_do))
        if should_save == 1
            rho_kr = cell(K_se, N_ops);
        end
        beta_kr = zeros(K_se, N_ops);
        beta_kr_time = zeros(K_se, N_ops);
    end
    
    %% Bosonic Symmetric Extensions
    %% Without PPT
    % Without realignment
    if any(ismember('beta_prime_konly', var_list))
        load(filename, 'beta_prime_konly', 'beta_prime_konly_time');
    elseif ~any(ismember('BosSymExtOnly', dont_do))
        if should_save == 1
            rho_prime_konly = cell(K_se, N_ops);
        end
        beta_prime_konly = zeros(K_se, N_ops);
        beta_prime_konly_time = zeros(K_se, N_ops);
    end

    % With realignment
    if any(ismember('beta_prime_konlyr', var_list))
        load(filename, 'beta_prime_konlyr', 'beta_prime_konlyr_time');
    elseif ~any(ismember('BostSymExtrlg', dont_do))
        if should_save == 1
            rho_prime_konlyr = cell(K_se, N_ops);
        end
        beta_prime_konlyr = zeros(K_se, N_ops);
        beta_prime_konlyr_time = zeros(K_se, N_ops);
    end

    %% With PPT
    % Without realignment
    if any(ismember('beta_prime_k', var_list))
        load(filename, 'beta_prime_k', 'beta_prime_k_time');
    elseif ~any(ismember('BosSymExt', dont_do))
        if should_save == 1
            rho_prime_k = cell(K_bse, N_ops);
        end
        beta_prime_k = zeros(K_se, N_ops);
        beta_prime_k_time = zeros(K_se, N_ops);
    end

    % With realignment
    if any(ismember('beta_prime_kr', var_list))
        load(filename, 'beta_prime_kr', 'beta_prime_kr_time');
    elseif ~any(ismember('BosSymExtrlg', dont_do))
        if should_save == 1
            rho_prime_kr = cell(K_bse, N_ops);
        end
        beta_prime_kr = zeros(K_se, N_ops);
        beta_prime_kr_time = zeros(K_se, N_ops);
    end
    
    %% Maximally mixed state
    % simple see-saw
    if any(ismember('gamma_MM', var_list))
        load(filename, 'gamma_MM', 'gamma_MM_stop', 'MM_time');
    elseif ~any(ismember('MM_simple', dont_do))
        if should_save == 1
            sigma_A_MM = cell(1, N_ops);
            sigma_B_MM = cell(1, N_ops);
        end
        gamma_MM = zeros(1, N_ops);
        gamma_MM_stop = zeros(1, N_ops);
        MM_time = zeros(1, N_ops);
    end
    
    % revised 1 see-saw
    if any(ismember('gamma_MM_rev_1', var_list))
        load(filename, 'gamma_MM_rev_1', 'gamma_MM_stop_rev_1', ...
                'MM_time_rev_1')
    elseif ~any(ismember('MM_rev_1', dont_do))
        if should_save == 1
            sigma_A_MM_rev_1 = cell(1, N_ops);
            sigma_B_MM_rev_1 = cell(1, N_ops);
        end
        gamma_MM_rev_1 = zeros(1, N_ops);
        gamma_MM_stop_rev_1 = zeros(1, N_ops);
        MM_time_rev_1 = zeros(1, N_ops);
    end
    
    % revised 2 see-saw
    if any(ismember('gamma_MM_rev_2', var_list))
        load(filename, 'gamma_MM_rev_2', 'gamma_MM_stop_rev_2', ...
                'MM_time_rev_2');
    elseif ~any(ismember('MM_rev_2', dont_do))
        if should_save == 1
            sigma_A_MM_rev_2 = cell(1, N_ops);
            sigma_B_MM_rev_2 = cell(1, N_ops);
        end
        gamma_MM_rev_2 = zeros(1, N_ops);
        gamma_MM_stop_rev_2 = zeros(1, N_ops);
        MM_time_rev_2 = zeros(1, N_ops);
    end
    
    %% Uniform superposition state
    % simple see-saw
    if any(ismember('gamma_US', var_list))
        load(filename, 'gamma_US', 'gamma_US_stop', 'US_time');
    elseif ~any(ismember('US_simple', dont_do))
        if should_save == 1
            sigma_A_US = cell(1, N_ops);
            sigma_B_US = cell(1, N_ops);
        end
        gamma_US = zeros(1, N_ops);
        gamma_US_stop = zeros(1, N_ops);
        US_time = zeros(1, N_ops);
    end
    
    % revised 1 see-saw
    if any(ismember('gamma_US_rev_1', var_list))
        load(filename, 'gamma_US_rev_1', 'gamma_US_stop_rev_1', ...
                'US_time_rev_1')
    elseif ~any(ismember('US_rev_1', dont_do))
        if should_save == 1
            sigma_A_US_rev_1 = cell(1, N_ops);
            sigma_B_US_rev_1 = cell(1, N_ops);
        end
        gamma_US_rev_1 = zeros(1, N_ops);
        gamma_US_stop_rev_1 = zeros(1, N_ops);
        US_time_rev_1 = zeros(1, N_ops);
    end
    
    % revised 2 see-saw
    if any(ismember('gamma_US_rev_2', var_list))
        load(filename, 'gamma_US_rev_2', 'gamma_US_stop_rev_2', ...
                'US_time_rev_2');
    elseif ~any(ismember('US_rev_2', dont_do))
        if should_save == 1
            sigma_A_US_rev_2 = cell(1, N_ops);
            sigma_B_US_rev_2 = cell(1, N_ops);
        end
        gamma_US_rev_2 = zeros(1, N_ops);
        gamma_US_stop_rev_2 = zeros(1, N_ops);
        US_time_rev_2 = zeros(1, N_ops);
    end
    
    %% Random
    % simple see-saw
    if any(ismember('gamma_rand', var_list))
        load(filename, 'gamma_rand_list', 'gamma_rand', 'gamma_rand_stop', ...
            'rand_time');
    elseif ~any(ismember('rand_simple', dont_do))
        if should_save == 1
            sigma_A_rand = cell(1, N_ops);
            sigma_B_rand = cell(1, N_ops);
        end
        gamma_rand_list = zeros(N_ops, N_rand+plus_two_rand);
        gamma_rand = zeros(1, N_ops);
        gamma_rand_stop = zeros(1, N_ops);
        rand_time = zeros(1, N_ops);
    end
    
    % revised 1 see-saw
    if any(ismember('gamma_rand_rev_1', var_list))
        load(filename, 'gamma_rand_list_rev_1', 'gamma_rand_rev_1', ...
            'gamma_rand_stop_rev_1', 'rand_time_rev_1');
    elseif ~any(ismember('rand_rev_1', dont_do))
        if should_save == 1
            sigma_A_rand_rev_1 = cell(1, N_ops);
            sigma_B_rand_rev_1 = cell(1, N_ops);
        end
        gamma_rand_list_rev_1 = zeros(N_ops, N_rand+plus_two_rand);
        gamma_rand_rev_1 = zeros(1, N_ops);
        gamma_rand_stop_rev_1 = zeros(1, N_ops);
        rand_time_rev_1 = zeros(1, N_ops);
    end
    
    % revised 2 see-saw
    if any(ismember('gamma_rand_rev_2', var_list))
        load(filename, 'gamma_rand_list_rev_2', 'gamma_rand_rev_2', ...
            'gamma_rand_stop_rev_2', 'rand_time_rev_2');
    elseif ~any(ismember('rand_rev_2', dont_do))
        if should_save == 1
            sigma_A_rand_rev_2 = cell(1, N_ops);
            sigma_B_rand_rev_2 = cell(1, N_ops);
        end
        gamma_rand_list_rev_2 = zeros(N_ops, N_rand+plus_two_rand);
        gamma_rand_rev_2 = zeros(1, N_ops);
        gamma_rand_stop_rev_2 = zeros(1, N_ops);
        rand_time_rev_2 = zeros(1, N_ops);
    end
        
    %% Tightness of bound
    ToB = zeros(1, N_ops);
    
    h = waitbar(0, "Computing approximations");
    for n = 1:N_ops
        %% Operator
        if sum == 1
            if ~any(ismember('K_list', var_list))
                genOtimesHouseholder(filename_ext, dim_A, dim_B, ...
                                     num_summands, opt);
            end
            load(filename, 'K_list', 'L_list');
%             Pi = zeros(dim_A * dim_B);
%             for i = 1:num_summands
%                 Pi = Pi + Tensor(K_list{i}, L_list{i});
%             end
        else
            if isempty(Pis{n})
                h = waitbar(0, "Generating operator");
                Pi = genMat(dim_A * dim_B, opt, proj);
                Pis{n} = Pi;
                %eig_Pis(:, n) = sort(eig(Pi), 'descend');
            else
                Pi = Pis{n};
            end
        end
        
        %% beta_PPT: PPT
        if ~any(ismember('PPT', dont_do))
            waitbar(n/N_ops, h, sprintf("Instance %d: PPT", n));
            PPT_start = tic;
            if sum == 1
                [rho_PPT_n, beta_PPT(n)] = computebeta_PPT_sum(K_list, ...
                                            L_list, num_summands);
            else
                [rho_PPT_n, beta_PPT(n)] = computebeta_PPT(Pi);
            end
            if should_save == 1
                rho_PPT{n} = rho_PPT_n;
            end
            beta_PPT_time(n) = toc(PPT_start);
            clear rho_PPT_n PPT_start;
        end
        
        %% beta_r: Realignment
        if ~any(ismember('realignment', dont_do))
            waitbar(n/N_ops, h, sprintf("Instance %d: Realignment", n));
            r_start = tic;
            if sum == 1
                [rho_r_n, beta_r(n)] = computebeta_r_sum(K_list, L_list, ...
                                                            num_summands);
            else
                [rho_r_n, beta_r(n)] = computebeta_r(Pi);
            end
            if should_save == 1
                rho_r{n} = rho_r_n;
            end
            beta_r_time(n) = toc(r_start);
            clear rho_r_n r_start;
        end
        
        %% beta_k: Symmetric Extensions
        %% Without PPT
        % Without realignment
        if ~any(ismember('SymExtOnly', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Symmetric Extensions (no PPT no realignment)", n));
            for k_se = 1:K_se
                konly_start = tic;
                if sum == 1
                    [rho_konly_se_n, ...
                        beta_konly(k_se, n)] = computebeta_k_sum(K_list, ...
                                                L_list, num_summands, ...
                                                dim_B, k_se, 0);
                else
                    [rho_konly_se_n, ...
                        beta_konly(k_se, n)] = computebeta_k(Pi, dim_B, ...
                                                                k_se, 0);
                end
                if should_save == 1
                    rho_konly{k_se}{n} = rho_konly_se_n;
                end
                beta_konly_time(k_se, n) = toc(konly_start);
                clear rho_konly_se_n konly_start;
            end
        end

        % With realignment
        if ~any(ismember('SymExtOnlyrlg', dont_do))
            waitbar(n/N_ops, h, ...
     sprintf("Instance %d: Symmetric Extensions realignment (no PPT)", n));
            for k_se = 1:K_se
                konlyr_start = tic;
                if sum == 1
                    [rho_konlyr_se_n, ...
                        beta_konlyr(k_se, n)] = computebeta_kr_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_B, k_se, 0);
                else
                    [rho_konlyr_se_n, ...
                        beta_konlyr(k_se, n)] = computebeta_kr(Pi, dim_B, ...
                                                                k_se, 0);
                end
                if should_save == 1
                    rho_konlyr{k_se}{n} = rho_konlyr_se_n;
                end
                beta_konlyr_time(k_se, n) = toc(konlyr_start);
                clear rho_konlyr_se_n konlyr_start;
            end
        end

        %% With PPT
        % Without realignment
        if ~any(ismember('SymExt', dont_do))
            waitbar(n/N_ops, h, ...
     sprintf("Instance %d: Symmetric Extensions PPT (no realignment)", n));
            for k_se = 1:K_se
                k_start = tic;
                if sum == 1
                    [rho_k_se_n, ...
                        beta_k(k_se, n)] = computebeta_k_sum(K_list, ...
                                            L_list, num_summands, dim_B, ...
                                            k_se, 1);
                else
                    [rho_k_se_n, beta_k(k_se, n)] = computebeta_k(Pi, ...
                                                        dim_B, k_se, 1);
                end
                if should_save == 1
                    rho_k{k_se}{n} = rho_k_se_n;
                end
                beta_k_time(k_se, n) = toc(k_start);
                clear rho_k_se_n k_start;
            end
        end

        % With realignment
        if ~any(ismember('SymExtrlg', dont_do))
            waitbar(n/N_ops, h, ...
          sprintf("Instance %d: Symmetric Extensions PPT realignment", n));
            for k_se = 1:K_se
                kr_start = tic;
                if sum == 1
                    [rho_kr_se_n, ...
                        beta_kr(k_se, n)] = computebeta_kr_sum(K_list, ...
                                                L_list, num_summands, ...
                                                dim_B, k_se, 1);
                else
                    [rho_kr_se_n, ...
                        beta_kr(k_se, n)] = computebeta_kr(Pi, dim_B, ...
                                                            k_se, 1);
                end
                if should_save == 1
                    rho_kr{k_se}{n} = rho_kr_se_n;
                end
                beta_kr_time(k_se, n) = toc(kr_start);
                clear rho_kr_se_n kr_start;
            end
        end
        
        %% beta_prime_k: Bosonic Symmetric Extensions
        %% Without PPT
        % Without realignment
        if ~any(ismember('BosSymExtOnly', dont_do))
            waitbar(n/N_ops, h, ...
sprintf("Instance %d: Bosonic Symmetric Extensions (no PPT no realignment)", ...
                    n));
            for k_bse = 1:K_bse
                prime_konly_start = tic;
                if sum == 1
                    [rho_prime_konly_bse_n, ...
                        beta_prime_konly(k_bse, n)] = computebeta_prime_k_sum(K_list, ...
                                                        L_list, num_summands, ...
                                                        dim_B, k_bse, 0);
                else
                    [rho_prime_konly_bse_n, ...
                        beta_prime_konly(k_bse, n)] = computebeta_prime_k(Pi, ...
                                                        dim_B, k_bse, 0);
                end
                if should_save == 1
                    rho_prime_konly{k_bse}{n} = rho_prime_konly_bse_n;
                end
                beta_prime_konly_time(k_bse, n) = toc(prime_konly_start);
                clear rho_prime_konly_bse_n prime_konly_start;
            end
        end

        % With realignment
        if ~any(ismember('BosSymExtOnlyrlg', dont_do))
            waitbar(n/N_ops, h, ...
sprintf("Instance %d: Bosonic Symmetric Extensions realignment (no PPT)", ...
                    n));
            for k_bse = 1:K_bse
                prime_konlyr_start = tic;
                if sum == 1
                    [rho_prime_konlyr_bse_n, ...
                        beta_prime_konlyr(k_bse, n)] = computebeta_prime_kr_sum(K_list, ...
                                                            L_list, ...
                                                            num_summands, ...
                                                            dim_B, k_bse, ...
                                                            0);
                else
                    [rho_prime_konlyr_bse_n, ...
                        beta_prime_konlyr(k_bse, n)] = computebeta_prime_kr(Pi, ...
                                                            dim_B, k_bse, 0);
                end
                if should_save == 1
                    rho_prime_konlyr{k_bse}{n} = rho_prime_konlyr_bse_n;
                end
                beta_prime_konlyr_time(k_bse, n) = toc(prime_konlyr_start);
                clear rho_prime_konlyr_bse_n prime_konlyr_start;
            end
        end

        %% With PPT
        % Without realignment
        if ~any(ismember('BosSymExt', dont_do))
            waitbar(n/N_ops, h, ...
sprintf("Instance %d: Bosonic Symmetric Extensions PPT (no realignment)", ...
                    n));
            for k_bse = 1:K_bse
                prime_k_start = tic;
                if sum == 1
                    [rho_prime_k_bse_n, ...
               beta_prime_k(k_bse, n)] = computebeta_prime_k_sum(K_list, ...
                                            L_list, num_summands, dim_B, ...
                                            k_bse, 1);
                else
                    [rho_prime_k_bse_n, ...
                        beta_prime_k(k_bse, n)] = computebeta_prime_k(Pi, ...
                                                            dim_B, k_bse, ...
                                                            1);
                end
                if should_save == 1
                    rho_prime_k{k_bse}{n} = rho_prime_k_bse_n;
                end
                beta_prime_k_time(k_bse, n) = toc(prime_k_start);
                clear rho_prime_k_bse_n prime_k_start;
            end
        end
        
        % With realignment
        if ~any(ismember('BosSymExtrlg', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Bosonic Symmetric Extensions PPT realignment", n));
            for k_bse = 1:K_bse
                prime_kr_start = tic;
                if sum == 1
                    [rho_prime_kr_bse_n, ...
                        beta_prime_kr(k_bse, n)] = computebeta_prime_kr_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_B, k_bse, 1);
                else
                    [rho_prime_kr_bse_n, ...
                        beta_prime_kr(k_bse, n)] = compute_prime_kr(Pi, ...
                                                    dim_B, k_bse, 1);
                end
                if should_save == 1
                    rho_prime_kr{k_bse}{n} = rho_prime_kr_bse_n;
                end
                beta_prime_kr_time(k_bse, n) = toc(prime_kr_start);
                clear rho_prime_kr_bse_n prime_kr_start;
            end
        end

        %% gamma_MM: See-saw from maximally-mixed state
        % simple see-saw
        if ~any(ismember('MM_simple', dont_do))
            waitbar(n/N_ops, h, ...
        sprintf("Instance %d: Maximally-mixed state (simple see-saw)", n));
            MM_start = tic;
            if sum == 1
                [sigma_A_MM_n, sigma_B_MM_n, gamma_MM(n), ...
                    gamma_MM_stop(n)] = computegamma_part_sum(K_list, ...
                                            L_list, num_summands, ...
                                            eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            else
                [sigma_A_MM_n, sigma_B_MM_n, gamma_MM(n), ...
                    gamma_MM_stop(n)] = computegamma_part(Pi, ...
                                                   eye(dim_A)/dim_A, ...
                                                   dim_B, threshold, ...
                                                   N_seeSaw, "sigma_A");
            end
            if should_save == 1
                sigma_A_MM{n} = sigma_A_MM_n;
                sigma_B_MM{n} = sigma_B_MM_n;
            end
            MM_time(n) = toc(MM_start);
            clear sigma_A_MM_n sigma_B_MM_n MM_start;
        end
        
        % revised 1 see-saw
        if ~any(ismember('MM_rev_1', dont_do))
            waitbar(n/N_ops, h, ...
     sprintf("Instance %d: Maximally-mixed state (revised 1 see-saw)", n));
            MM_start_rev_1 = tic;
            if sum == 1
                [sigma_A_MM_rev_1_n, sigma_B_MM_rev_1_n, ...
                    gamma_MM_rev_1(n), ...
                 gamma_MM_stop_rev_1(n)] = computegamma_part_rev_1_sum(K_list, ...
                                            L_list, num_summands, ...
                                            eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            else
                [sigma_A_MM_rev_1_n, sigma_B_MM_rev_1_n, ...
                    gamma_MM_rev_1(n), ...
                 gamma_MM_stop_rev_1(n)] = computegamma_part_rev_1(Pi, ...
                                            eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            end
            if should_save == 1
                sigma_A_MM_rev_1{n} = sigma_A_MM_rev_1_n;
                sigma_B_MM_rev_1{n} = sigma_B_MM_rev_1_n;
            end
            MM_time_rev_1(n) = toc(MM_start_rev_1);
            clear sigma_A_MM_rev_1_n sigma_B_MM_rev_1_n MM_start_rev_1;
        end
        
        % revised 2 see-saw
        if ~any(ismember('MM_rev_2', dont_do))
            waitbar(n/N_ops, h, ...
     sprintf("Instance %d: Maximally-mixed state (revised 2 see-saw)", n));
            MM_start_rev_2 = tic;
            if sum == 1
                [sigma_A_MM_rev_2_n, sigma_B_MM_rev_2_n, ...
                    gamma_MM_rev_2(n), ...
                 gamma_MM_stop_rev_2(n)] = computegamma_part_rev_2_sum(K_list, ...
                                            L_list, num_summands, ...
                                            eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            else
                [sigma_A_MM_rev_2_n, sigma_B_MM_rev_2_n, ...
                    gamma_MM_rev_2(n), ...
                 gamma_MM_stop_rev_2(n)] = computegamma_part_rev_2(Pi, ...
                                            eye(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            end
            if should_save == 1
                sigma_A_MM_rev_2{n} = sigma_A_MM_rev_2_n;
                sigma_B_MM_rev_2{n} = sigma_B_MM_rev_2_n;
            end
            MM_time_rev_2(n) = toc(MM_start_rev_2);
            clear sigma_A_MM_rev_2_n sigma_B_MM_rev_2_n MM_start_rev_2;
        end
        
        %% gamma_US: See-saw from uniform superposition state
        % simple see-saw
        if ~any(ismember('US_simple', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Uniform superposition state (simple see-saw)", n));
            US_start = tic;
            if sum == 1
                [sigma_A_US_n, sigma_B_US_n, gamma_US(n), ...
                    gamma_US_stop(n)] = computegamma_part_sum(K_list, ...
                                            L_list, num_summands, ...
                                            ones(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            else
                [sigma_A_US_n, sigma_B_US_n, gamma_US(n), ...
                    gamma_US_stop(n)] = computegamma_part(Pi, ...
                                                   ones(dim_A)/dim_A, ...
                                                   dim_B, threshold, ...
                                                   N_seeSaw, "sigma_A");
            end
            if should_save == 1
                sigma_A_US{n} = sigma_A_US_n;
                sigma_B_US{n} = sigma_B_US_n;
            end
            US_time(n) = toc(US_start);
            clear sigma_A_US_n sigma_B_US_n US_start;
        end
        
        % revised 1 see-saw
        if ~any(ismember('US_rev_1', dont_do))
            waitbar(n/N_ops, h, ...
sprintf("Instance %d: Uniform superposition state (revised 1 see-saw)", n));
            US_start_rev_1 = tic;
            if sum == 1
                [sigma_A_US_rev_1_n, sigma_B_US_rev_1_n, ...
                    gamma_US_rev_1(n), ...
                    gamma_US_stop_rev_1(n)] = computegamma_part_rev_1_sum(K_list, ...
                                                L_list, num_summands, ...
                                                ones(dim_A)/dim_A, ...
                                                dim_B, threshold, ...
                                                N_seeSaw, "sigma_A");
            else
                [sigma_A_US_rev_1_n, sigma_B_US_rev_1_n, ...
                    gamma_US_rev_1(n), ...
                 gamma_US_stop_rev_1(n)] = computegamma_part_rev_1(Pi, ...
                                            ones(dim_A)/dim_A, dim_B, ...
                                            threshold, N_seeSaw, "sigma_A");
            end
            if should_save == 1
                sigma_A_US_rev_1{n} = sigma_A_US_rev_1_n;
                sigma_B_US_rev_1{n} = sigma_B_US_rev_1_n;
            end
            US_time_rev_1(n) = toc(US_start_rev_1);
            clear sigma_A_US_rev_1_n sigma_B_US_rev_1_n US_start_rev_1;
        end
        
        % revised 2 see-saw
        if ~any(ismember('US_rev_2', dont_do))
            waitbar(n/N_ops, h, ...
sprintf("Instance %d: Uniform superposition state (revised 2 see-saw)", n));
            US_start_rev_2 = tic;
            if sum == 1
                [sigma_A_US_rev_2_n, sigma_B_US_rev_2_n, ...
                    gamma_US_rev_2(n), ...
                    gamma_US_stop_rev_2(n)] = computegamma_part_rev_2_sum(K_list, ...
                                                L_list, num_summands, ...
                                                ones(dim_A)/dim_A, ...
                                                dim_B, threshold, ...
                                                N_seeSaw, "sigma_A");
            else
                [sigma_A_US_rev_2_n, sigma_B_US_rev_2_n, ...
                    gamma_US_rev_2(n), ...
                    gamma_US_stop_rev_2(n)] = computegamma_part_rev_2(Pi, ...
                                                ones(dim_A)/dim_A, dim_B, ...
                                                threshold, N_seeSaw, ...
                                                "sigma_A");
            end
            if should_save == 1
                sigma_A_US_rev_2{n} = sigma_A_US_rev_2_n;
                sigma_B_US_rev_2{n} = sigma_B_US_rev_2_n;
            end
            US_time_rev_2(n) = toc(US_start_rev_2);
            clear sigma_A_US_rev_2_n sigma_B_US_rev_2_n US_start_rev_2;
        end
        
        %% gamma_rand: Random see-saw
        % simple see-saw
        if ~any(ismember('rand_simple', dont_do))
            waitbar(n/N_ops, h, ...
                    sprintf("Instance %d: Random (simple see-saw)", n));
            rand_start = tic;
            if sum == 1
                [sigma_A_rand_n, sigma_B_rand_n, gamma_rand(n), ...
                    gamma_rand_stop(n), ...
                    gamma_rand_list(n,:)] = computegamma_rand_sum(K_list, ...
                                                L_list, num_summands, dim_A, ...
                                                dim_B, threshold, N_seeSaw, ...
                                                N_rand, "simple", ...
                                                plus_two_rand);
            else
                [sigma_A_rand_n, sigma_B_rand_n, gamma_rand(n), ...
                    gamma_rand_stop(n), ...
                    gamma_rand_list(n,:)] = computegamma_rand(Pi, dim_A, ...
                                                dim_B, threshold, N_seeSaw, ...
                                                N_rand, "simple", ...
                                                plus_two_rand);
            end
            if should_save == 1
                sigma_A_rand{n} = sigma_A_rand_n;
                sigma_B_rand{n} = sigma_B_rand_n;
            end
            rand_time(n) = toc(rand_start);
            clear sigma_A_rand_n sigma_B_rand_n rand_start;
        end
        
        % revised 1 see-saw
        if ~any(ismember('rand_rev_1', dont_do))
            waitbar(n/N_ops, h, ...
                    sprintf("Instance %d: Random (revised 1 see-saw)", n));
            rand_start_rev_1 = tic;
            if sum == 1
                [sigma_A_rand_rev_1_n, sigma_B_rand_rev_1_n, ...
                    gamma_rand_rev_1(n), gamma_rand_stop_rev_1(n), ...
                    gamma_rand_list_rev_1(n,:)] = computegamma_rand_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_A, dim_B, threshold, ...
                                                    N_seeSaw, N_rand, ...
                                                    "rev_1", plus_two_rand);
            else
                [sigma_A_rand_rev_1_n, sigma_B_rand_rev_1_n, ...
                    gamma_rand_rev_1(n), gamma_rand_stop_rev_1(n), ...
                    gamma_rand_list_rev_1(n,:)] = computegamma_rand(Pi, ...
                                                    dim_A, dim_B, threshold, ...
                                                    N_seeSaw, N_rand, ...
                                                    "rev_1", plus_two_rand);
            end
            if should_save == 1
                sigma_A_rand_rev_1{n} = sigma_A_rand_rev_1_n;
                sigma_B_rand_rev_1{n} = sigma_B_rand_rev_1_n;
            end
            rand_time_rev_1(n) = toc(rand_start_rev_1);
            clear sigma_A_rand_rev_1_n sigma_B_rand_rev_1_n ...
                rand_start_rev_1;
        end
        
        % revised 2 see-saw
        if ~any(ismember('rand_rev_2', dont_do))
            waitbar(n/N_ops, h, ...
                    sprintf("Instance %d: Random (revised 2 see-saw)", n));
            rand_start_rev_2 = tic;
            if sum == 1
                [sigma_A_rand_rev_2_n, sigma_B_rand_rev_2_n, ...
                    gamma_rand_rev_2(n), gamma_rand_stop_rev_2(n), ...
                    gamma_rand_list_rev_2(n,:)] = computegamma_rand_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_A, dim_B, threshold, ...
                                                    N_seeSaw, N_rand, ...
                                                    "rev_2", plus_two_rand);
            else
                [sigma_A_rand_rev_2_n, sigma_B_rand_rev_2_n, ...
                    gamma_rand_rev_2(n), gamma_rand_stop_rev_2(n), ...
                    gamma_rand_list_rev_2(n,:)] = computegamma_rand(Pi, ...
                                                    dim_A, dim_B, threshold, ...
                                                    N_seeSaw, N_rand, ...
                                                    "rev_2", plus_two_rand);
            end
            if should_save == 1
                sigma_A_rand_rev_2{n} = sigma_A_rand_rev_2_n;
                sigma_B_rand_rev_2{n} = sigma_B_rand_rev_2_n;
            end
            rand_time_rev_2(n) = toc(rand_start_rev_2);
            clear sigma_A_rand_rev_2_n sigma_B_rand_rev_2_n ...
                rand_start_rev_2;
        end
        
        %% beta: Smallest upper bound
        beta = nan;
        
        % PPT
        if ~any(ismember('PPT', dont_do)) || any(ismember('beta_PPT', ...
                var_list))
            beta = min(beta, beta_PPT(n));
        end

        % Realignment
        if ~any(ismember('realignment', ...
                dont_do)) || any(ismember('beta_r', var_list))
            beta = min(beta, beta_r(n));
        end
        
        % Symmetric Extensions
        if ~any(ismember('SymExtOnly', ...
                dont_do)) || any(ismember('beta_konly', var_list))
            beta = min(beta, min(beta_konly(:, n)));
        end
        if ~any(ismember('SymExtOnlyrlg', ...
                dont_do)) || any(ismember('beta_konlyr', var_list))
            beta = min(beta, min(beta_konlyr(:, n)));
        end
        if ~any(ismember('SymExt', ...
                dont_do)) || any(ismember('beta_k', var_list))
            beta = min(beta, min(beta_k(:, n)));
        end
        if ~any(ismember('SymExtrlg', ...
                dont_do)) || any(ismember('beta_kr', var_list))
            beta = min(beta, min(beta_kr(:, n)));
        end
        
        % Bosonic Symmetric Extensions
        if ~any(ismember('BosSymExtOnly', ...
                dont_do)) || any(ismember('beta_prime_konly', var_list))
            beta = min(beta, min(beta_prime_konly(:, n)));
        end
        if ~any(ismember('BosSymExtOnlyrlg', ...
                dont_do)) || any(ismember('beta_prime_konlyr', var_list))
            beta = min(beta, min(beta_prime_konlyr(:, n)));
        end
        if ~any(ismember('BosSymExt', ...
                dont_do)) || any(ismember('beta_prime_k', var_list))
            beta = min(beta, min(beta_prime_k(:, n)));
        end
        if ~any(ismember('BosSymExtrlg', ...
                dont_do)) || any(ismember('beta_prime_kr', var_list))
            beta = min(beta, min(beta_prime_kr(:, n)));
        end
        
        %% gamma: Largest upper bound
        gamma = nan;
        
        % Maximally-mixed state
        if ~any(ismember('MM_simple', ...
                dont_do)) || any(ismember('gamma_MM', var_list))
            gamma = max(gamma, gamma_MM(n));
        end
        if ~any(ismember('MM_rev_1', ...
                dont_do)) || any(ismember('gamma_MM_rev_1', var_list))
            gamma = max(gamma, gamma_MM_rev_1(n));
        end
        if ~any(ismember('MM_rev_2', ...
                dont_do)) || any(ismember('gamma_MM_rev_2', var_list))
            gamma = max(gamma, gamma_MM_rev_2(n));
        end
        
        % Uniform superposition state
        if ~any(ismember('US_simple', ...
                dont_do)) || any(ismember('gamma_US', var_list))
            gamma = max(gamma, gamma_US(n));
        end
        if ~any(ismember('US_rev_1', ...
                dont_do)) || any(ismember('gamma_US_rev_1', var_list))
            gamma = max(gamma, gamma_US_rev_1(n));
        end
        if ~any(ismember('US_rev_2', ...
                dont_do)) || any(ismember('gamma_US_rev_2', var_list))
            gamma = max(gamma, gamma_US_rev_2(n));
        end
        
        % Random starting points
        if ~any(ismember('rand_simple', ...
                dont_do)) || any(ismember('gamma_rand', var_list))
            gamma = max(gamma, gamma_rand(n));
        end
        if ~any(ismember('rand_rev_1', ...
                dont_do)) || any(ismember('gamma_rand_rev_1', var_list))
            gamma = max(gamma, gamma_rand_rev_1(n));
        end
        if ~any(ismember('rand_rev_2', ...
                dont_do)) || any(ismember('gamma_rand_rev_2', var_list))
            gamma = max(gamma, gamma_rand_rev_2(n));
        end
        
        %% Tightness of bound
        waitbar(n/N_ops, h, sprintf("Instance %d: Computing ToB", n));
        ToB(n) = beta - gamma;
    end
    %% Save arrays to filename
    waitbar(n/N_ops, h, sprintf("Instance %d: Saving arrays", n));
    
    if ~isfile(filename_ext)
        if sum == 0
            save(filename, 'Pis', '-v7.3');
            clear Pis;
        end
    end

    % PPT
    if ~any(ismember('PPT', dont_do))
        if should_save == 1
            save(filename, 'rho_PPT', '-append', '-v7.3');
            clear rho_PPT;
        end
        save(filename, 'beta_PPT', 'beta_PPT_time', '-append', '-v7.3');
        clear beta_PPT beta_PPT_time;
    end

    % Realignment
    if ~any(ismember('realignment', dont_do))
        if should_save == 1
            save(filename, 'rho_r', '-append', '-v7.3');
            clear rho_r;
        end
        save(filename, 'beta_r', 'beta_r_time', '-append', '-v7.3');
        clear beta_r beta_r_time;
    end

    % Symmetric Extensions
    if ~any(ismember('SymExtOnly', dont_do))
        if should_save == 1
            save(filename, 'rho_konly', '-append', '-v7.3');
            clear rho_konly;
        end
        save(filename, 'beta_konly', 'beta_konly_time', '-append', ...
            '-v7.3');
        clear beta_konly beta_konly_time;
    end
    if ~any(ismember('SymExtOnlyrlg', dont_do))
        if should_save == 1
            save(filename, 'rho_konlyr', '-append', '-v7.3');
            clear rho_konlyr;
        end
        save(filename, 'beta_konlyr', 'beta_konlyr_time', '-append', ...
            '-v7.3');
        clear beta_konlyr beta_konlyr_time;
    end
    if ~any(ismember('SymExt', dont_do))
        if should_save == 1
            save(filename, 'rho_k', '-append', '-v7.3');
            clear rho_k;
        end
        save(filename, 'beta_k', 'beta_k_time', '-append', '-v7.3');
        clear beta_k beta_k_time;
    end
    if ~any(ismember('SymExtrlg', dont_do))
        if should_save == 1
            save(filename, 'rho_kr', '-append', '-v7.3');
            clear rho_kr;
        end
        save(filename, 'beta_kr', 'beta_kr_time', '-append', '-v7.3');
        clear beta_kr beta_kr_time;
    end

    % Bosonic Symmetric Extensions
    if ~any(ismember('BosSymExtOnly', dont_do))
        if should_save == 1
            save(filename, 'rho_prime_konly', '-append', '-v7.3');
            clear rho_prime_konly;
        end
        save(filename, 'beta_prime_konly', 'beta_prime_konly_time', ...
            '-append', '-v7.3');
        clear beta_prime_konly beta_prime_konly_time;
    end
    if ~any(ismember('BosSymExtOnlyrlg', dont_do))
        if should_save == 1
            save(filename, 'rho_prime_konlyr', '-append', '-v7.3');
            clear rho_prime_konlyr;
        end
        save(filename, 'beta_prime_konlyr', 'beta_prime_konlyr_time', ...
            '-append', '-v7.3');
        clear beta_prime_konlyr beta_prime_konlyr_time;
    end
    if ~any(ismember('BosSymExt', dont_do))
        if should_save == 1
            save(filename, 'rho_prime_k', '-append', '-v7.3');
            clear rho_prime_k;
        end
        save(filename, 'beta_prime_k', 'beta_prime_k_time', '-append', ...
                '-v7.3');
        clear beta_prime_k beta_prime_k_time;
    end
    if ~any(ismember('BosSymExtrlg', dont_do))
        if should_save == 1
            save(filename, 'rho_prime_kr', '-append', '-v7.3');
            clear rho_prime_kr;
        end
        save(filename, 'beta_prime_kr', 'beta_prime_kr_time', '-append', ...
                '-v7.3');
        clear beta_prime_kr beta_prime_kr_time;
    end
    
    % Maximally-mixed state
    if ~any(ismember('MM_simple', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_MM', 'sigma_B_MM', '-append', '-v7.3');
            clear sigma_A_MM sigma_B_MM;
        end
        save(filename, 'gamma_MM', 'gamma_MM_stop', 'MM_time', ...
                '-append', '-v7.3');
        clear gamma_MM gamma_MM_stop MM_time;
    end
    if ~any(ismember('MM_rev_1', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_MM_rev_1', 'sigma_B_MM_rev_1', ...
                    '-append', '-v7.3');
            clear sigma_A_MM_rev_1 sigma_B_MM_rev_1;
        end
        save(filename, 'gamma_MM_rev_1', 'gamma_MM_stop_rev_1', ...
                'MM_time_rev_1', '-append', '-v7.3');
        clear gamma_MM_rev_1 gamma_MM_stop_rev_1 MM_time_rev_1;
    end
    if ~any(ismember('MM_rev_2', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_MM_rev_2', 'sigma_B_MM_rev_2', ...
                    '-append', '-v7.3');
            clear sigma_A_MM_rev_2 sigma_B_MM_rev_2;
        end
        save(filename, 'gamma_MM_rev_2', 'gamma_MM_stop_rev_2', ...
                'MM_time_rev_2', '-append', '-v7.3');
        clear gamma_MM_rev_2 gamma_MM_stop_rev_2 MM_time_rev_2;
    end
    
    % Uniform superposition state
    if ~any(ismember('US_simple', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_US', 'sigma_B_US', '-append', '-v7.3');
            clear sigma_A_US sigma_B_US;
        end
        save(filename, 'gamma_US', 'gamma_US_stop', 'US_time', ...
                '-append', '-v7.3');
        clear gamma_US gamma_US_stop US_time;
    end
    if ~any(ismember('US_rev_1', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_US_rev_1', 'sigma_B_US_rev_1', ...
                    '-append', '-v7.3');
            clear sigma_A_US_rev_1 sigma_B_US_rev_1;
        end
        save(filename, 'gamma_US_rev_1', 'gamma_US_stop_rev_1', ...
                'US_time_rev_1', '-append', '-v7.3');
        clear gamma_US_rev_1 gamma_US_stop_rev_1 US_time_rev_1;
    end
    if ~any(ismember('US_rev_2', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_US_rev_2', 'sigma_B_US_rev_2', ...
                    '-append', '-v7.3');
            clear sigma_A_US_rev_2 sigma_B_US_rev_2;
        end
        save(filename, 'gamma_US_rev_2', 'gamma_US_stop_rev_2', ...
                'US_time_rev_2', '-append', '-v7.3');
        clear gamma_US_rev_2 gamma_US_stop_rev_2 US_time_rev_2;
    end
    
    % Random starting point
    if ~any(ismember('rand_simple', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rand', 'sigma_B_rand', '-append', ...
                '-v7.3');
            clear sigma_A_rand sigma_B_rand;
        end
        save(filename, 'gamma_rand_list', 'gamma_rand', 'gamma_rand_stop', ...
            'rand_time', '-append', '-v7.3');
        clear gamma_rand_list gamma_rand gamma_rand_stop rand_time;
    end
    if ~any(ismember('rand_rev_1', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rand_rev_1', 'sigma_B_rand_rev_1', ...
                    '-append', '-v7.3');
            clear sigma_A_rand_rev_1 sigma_B_rand_rev_1;
        end
        save(filename, 'gamma_rand_list_rev_1', 'gamma_rand_rev_1', ...
            'gamma_rand_stop_rev_1', 'rand_time_rev_1', '-append', '-v7.3');
        clear gamma_rand_list_rev_1 gamma_rand_rev_1 ...
            gamma_rand_stop_rev_1 rand_time_rev_1;
    end
    if ~any(ismember('rand_rev_2', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rand_rev_2', 'sigma_B_rand_rev_2', ...
                    '-append', '-v7.3');
            clear sigma_A_rand_rev_2 sigma_B_rand_rev_2;
        end
        save(filename, 'gamma_rand_list_rev_2', 'gamma_rand_rev_2', ...
            'gamma_rand_stop_rev_2', 'rand_time_rev_2', '-append', '-v7.3');
        clear gamma_rand_list_rev_2 gamma_rand_rev_2 ...
            gamma_rand_stop_rev_2 rand_time_rev_2;
    end
    
    save(filename, 'ToB', '-append', '-v7.3');
    clear ToB;
    
    close(h);
end