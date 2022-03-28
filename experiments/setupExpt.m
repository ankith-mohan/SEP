function setupExpt(dim_A, dim_B, K_se, K_bse, N_seeSaw, N_rand, ...
                   sum, num_summands, N_ops, threshold, opt, filename, ...
                   proj, dont_do, should_save)
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
    if any(ismember('beta_k', var_list))
        load(filename, 'beta_k', 'beta_k_time');
    elseif ~any(ismember('SymExt', dont_do))
        if should_save == 1
            rho_k = cell(K_se, N_ops);
        end
        beta_k = zeros(K_se, N_ops);
        beta_k_time = zeros(K_se, N_ops);
    end
    
    %% Bosonic Symmetric Extensions
    if any(ismember('beta_prime_k', var_list))
        load(filename, 'beta_prime_k', 'beta_prime_k_time');
    elseif ~any(ismember('BosSymExt', dont_do))
        if should_save == 1
            rho_prime_k = cell(K_bse, N_ops);
        end
        beta_prime_k = zeros(K_se, N_ops);
        beta_prime_k_time = zeros(K_se, N_ops);
    end
    
    %% Square root
    % simple see-saw from sigma_A
    if any(ismember('beta_sqrt', var_list))
        load(filename, 'beta_sqrt', 'gamma_sqrt', 'gamma_sqrt_stop', ...
                'sqrt_time');
    elseif ~any(ismember('sqrt_simple_A', dont_do))
        if should_save == 1
            sigma_A = cell(1, N_ops);
            sigma_B = cell(1, N_ops);
        end
        beta_sqrt = zeros(1, N_ops);
        
        if should_save == 1
            sigma_A_sqrt = cell(1, N_ops);
            sigma_B_sqrt = cell(1, N_ops);
        end
        gamma_sqrt = zeros(1, N_ops);
        
        gamma_sqrt_stop = zeros(1, N_ops);
        sqrt_time = zeros(1, N_ops);
    end
    
    % simple see-saw from sigma_B
    if any(ismember('beta_sqrt_B', var_list))
        load(filename, 'beta_sqrt_B', 'gamma_sqrt_B', ...
                'gamma_sqrt_stop_B', 'sqrt_time_B');
    elseif ~any(ismember('sqrt_simple_B', dont_do))
        if should_save == 1
            sigma_A_start_B = cell(1, N_ops);
            sigma_B_start_B = cell(1, N_ops);
        end
        beta_sqrt_B = zeros(1, N_ops);
        
        if should_save == 1
            sigma_A_sqrt_B = cell(1, N_ops);
            sigma_B_sqrt_B = cell(1, N_ops);
        end
        gamma_sqrt_B = zeros(1, N_ops);
        
        gamma_sqrt_stop_B = zeros(1, N_ops);
        sqrt_time_B = zeros(1, N_ops);
    end
    
    % revised 1 see-saw from sigma_A
    if any(ismember('beta_sqrt_rev_1', var_list))
        load(filename, 'beta_sqrt_rev_1', 'gamma_sqrt_rev_1', ...
                'gamma_sqrt_stop_rev_1', 'sqrt_time_rev_1');
    elseif ~any(ismember('sqrt_rev_1_A', dont_do))
        if should_save == 1
            sigma_A_rev_1 = cell(1, N_ops);
            sigma_B_rev_1 = cell(1, N_ops);
        end
        beta_sqrt_rev_1 = zeros(1, N_ops);
        
        if should_save == 1
            sigma_A_sqrt_rev_1 = cell(1, N_ops);
            sigma_B_sqrt_rev_1 = cell(1, N_ops);
        end
        gamma_sqrt_rev_1 = zeros(1, N_ops);
        
        gamma_sqrt_stop_rev_1 = zeros(1, N_ops);
        sqrt_time_rev_1 = zeros(1, N_ops);
    end
    
    % revised 1 see-saw from sigma_B
    if any(ismember('beta_sqrt_rev_1_B', var_list))
        load(filename, 'beta_sqrt_rev_1_B', 'gamma_sqrt_rev_1_B', ...
                'gamma_sqrt_stop_rev_1_B', 'sqrt_time_rev_1_B');
    elseif ~any(ismember('sqrt_rev_1_B', dont_do))
        if should_save == 1
            sigma_A_rev_1_B = cell(1, N_ops);
            sigma_B_rev_1_B = cell(1, N_ops);
        end
        beta_sqrt_rev_1_B = zeros(1, N_ops);
        
        if should_save == 1
            sigma_A_sqrt_rev_1_B = cell(1, N_ops);
            sigma_B_sqrt_rev_1_B = cell(1, N_ops);
        end
        gamma_sqrt_rev_1_B = zeros(1, N_ops);
        
        gamma_sqrt_stop_rev_1_B = zeros(1, N_ops);
        sqrt_time_rev_1_B = zeros(1, N_ops);
    end
    
    % revised 2 see-saw from sigma_A
    if any(ismember('beta_sqrt_rev_2', var_list))
        load(filename, 'beta_sqrt_rev_2', 'gamma_sqrt_rev_2', ...
                'gamma_sqrt_stop_rev_2', 'sqrt_time_rev_2');
    elseif ~any(ismember('sqrt_rev_2_A', dont_do))
        if should_save == 1
            sigma_A_rev_2 = cell(1, N_ops);
            sigma_B_rev_2 = cell(1, N_ops);
        end
        beta_sqrt_rev_2 = zeros(1, N_ops);
        
        if should_save == 1
            sigma_A_sqrt_rev_2 = cell(1, N_ops);
            sigma_B_sqrt_rev_2 = cell(1, N_ops);
        end
        gamma_sqrt_rev_2 = zeros(1, N_ops);
        
        gamma_sqrt_stop_rev_2 = zeros(1, N_ops);
        sqrt_time_rev_2 = zeros(1, N_ops);
    end
    
    % revised 2 see-saw from sigma_B
    if any(ismember('beta_sqrt_rev_2_B', var_list))
        load(filename, 'beta_sqrt_rev_2_B', 'gamma_sqrt_rev_2_B', ...
                'gamma_sqrt_stop_rev_2_B', 'sqrt_time_rev_2_B');
    elseif ~any(ismember('sqrt_rev_2_B', dont_do))
        if should_save == 1
            sigma_A_rev_2_B = cell(1, N_ops);
            sigma_B_rev_2_B = cell(1, N_ops);
        end
        beta_sqrt_rev_2_B = zeros(1, N_ops);
        
        if should_save == 1
            sigma_A_sqrt_rev_2_B = cell(1, N_ops);
            sigma_B_sqrt_rev_2_B = cell(1, N_ops);
        end
        gamma_sqrt_rev_2_B = zeros(1, N_ops);
        
        gamma_sqrt_stop_rev_2_B = zeros(1, N_ops);
        sqrt_time_rev_2_B = zeros(1, N_ops);
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
        load(filename, 'gamma_rand', 'gamma_rand_stop', 'rand_time');
    elseif ~any(ismember('rand_simple', dont_do))
        if should_save == 1
            sigma_A_rand = cell(1, N_ops);
            sigma_B_rand = cell(1, N_ops);
        end
        gamma_rand = zeros(1, N_ops);
        gamma_rand_stop = zeros(1, N_ops);
        rand_time = zeros(1, N_ops);
    end
    
    % revised 1 see-saw
    if any(ismember('gamma_rand_rev_1', var_list))
        load(filename, 'gamma_rand_rev_1', 'gamma_rand_stop_rev_1', ...
                'rand_time_rev_1')
    elseif ~any(ismember('rand_rev_1', dont_do))
        if should_save == 1
            sigma_A_rand_rev_1 = cell(1, N_ops);
            sigma_B_rand_rev_1 = cell(1, N_ops);
        end
        gamma_rand_rev_1 = zeros(1, N_ops);
        gamma_rand_stop_rev_1 = zeros(1, N_ops);
        rand_time_rev_1 = zeros(1, N_ops);
    end
    
    % revised 2 see-saw
    if any(ismember('gamma_rand_rev_2', var_list))
        load(filename, 'gamma_rand_rev_2', 'gamma_rand_stop_rev_2', ...
                'rand_time_rev_2');
    elseif ~any(ismember('rand_rev_2', dont_do))
        if should_save == 1
            sigma_A_rand_rev_2 = cell(1, N_ops);
            sigma_B_rand_rev_2 = cell(1, N_ops);
        end
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
        end
        
        %% beta_k: Symmetric Extensions
        if ~any(ismember('SymExt', dont_do))
            waitbar(n/N_ops, h, ...
                sprintf("Instance %d: Symmetric Extensions", n));
            for k_se = 1:K_se
                k_start = tic;
                if sum == 1
                    [rho_k_se_n, ...
                        beta_k(k_se, n)] = computebeta_k_sum(K_list, ...
                                            L_list, num_summands, dim_B, ...
                                            k_se);
                else
                    [rho_k_se_n, beta_k(k_se, n)] = computebeta_k(Pi, ...
                                                        dim_B, k_se);
                end
                if should_save == 1
                    rho_k{k_se}{n} = rho_k_se_n;
                end
                beta_k_time(k_se, n) = toc(k_start);
            end
        end
        
        %% beta_prime_k: Bosonic Symmetric Extensions
        if ~any(ismember('BosSymExt', dont_do))
            waitbar(n/N_ops, h, ...
                sprintf("Instance %d: Bosonic Symmetric Extensions", n));
            for k_bse = 1:K_bse
                prime_k_start = tic;
                if sum == 1
                    [rho_prime_k_bse_n, ...
               beta_prime_k(k_bse, n)] = computebeta_prime_k_sum(K_list, ...
                                            L_list, num_summands, dim_B, ...
                                            k_bse);
                else
                    [rho_prime_k_bse_n, ...
                        beta_prime_k(k_bse, n)] = computebeta_prime_k(Pi, ...
                                                            dim_B, k_bse);
                end
                if should_save == 1
                    rho_prime_k{k_bse}{n} = rho_prime_k_bse_n;
                end
                beta_prime_k_time(k_bse, n) = toc(prime_k_start);
            end
        end
        
        %% sqrt: Square root
        % beta_sqrt: Upper bound
        % gamma_sqrt: Lower bound
        
        % simple see-saw from sigma_A
        if ~any(ismember('sqrt_simple_A', dont_do))
            waitbar(n/N_ops, h, ...
     sprintf("Instance %d: Square root (simple see-saw from sigma_A)", n));
            sqrt_start = tic;
            if sum == 1
                [sigma_A_n, sigma_B_n, beta_sqrt(n), sigma_A_sqrt_n, ...
                    sigma_B_sqrt_n, gamma_sqrt(n), ...
                    gamma_sqrt_stop(n)] = computesqrt_sum(K_list, L_list, ...
                                        num_summands, dim_A, dim_B, ...
                                        threshold, N_seeSaw, "sigma_A", ...
                                        "simple");
            else
            [sigma_A_n, sigma_B_n, beta_sqrt(n), sigma_A_sqrt_n, ...
                sigma_B_sqrt_n, gamma_sqrt(n), ...
                gamma_sqrt_stop(n)] = computesqrt(Pi, dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    "sigma_A", "simple");
            end
            if should_save == 1
                sigma_A{n} = sigma_A_n;
                sigma_B{n} = sigma_B_n;
                sigma_A_sqrt{n} = sigma_A_sqrt_n;
                sigma_B_sqrt{n} = sigma_B_sqrt_n;
            end
            sqrt_time(n) = toc(sqrt_start);
        end
        
        % simple see-saw from sigma_B
        if ~any(ismember('sqrt_simple_B', dont_do))
            waitbar(n/N_ops, h, ...
     sprintf("Instance %d: Square root (simple see-saw from sigma_B)", n));
            sqrt_start_B = tic;
            if sum == 1
                [sigma_A_start_B_n, sigma_B_start_B_n, beta_sqrt_B(n), ...
                    sigma_A_sqrt_B_n, sigma_B_sqrt_B_n, gamma_sqrt_B(n), ...
                    gamma_sqrt_stop_B(n)] = computesqrt_sum(K_list, ...
                                                L_list, num_summands, ...
                                                dim_A, dim_B, threshold, ...
                                                N_seeSaw, "sigma_B", ...
                                                "simple"); 
            else
                [sigma_A_start_B_n, sigma_B_start_B_n, beta_sqrt_B(n), ...
                    sigma_A_sqrt_B_n, sigma_B_sqrt_B_n, gamma_sqrt_B(n), ...
                    gamma_sqrt_stop_B(n)] = computesqrt(Pi, dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    "sigma_B", "simple"); 
            end
            if should_save == 1
                sigma_A_start_B{n} = sigma_A_start_B_n;
                sigma_B_start_B{n} = sigma_B_start_B_n;
                sigma_A_sqrt_B{n} = sigma_A_sqrt_B_n;
                sigma_B_sqrt_B{n} = sigma_B_sqrt_B_n;
            end
            sqrt_time_B(n) = toc(sqrt_start_B);
        end
        
        % revised 1 see-saw from sigma_A
        if ~any(ismember('sqrt_rev_1_A', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Square root (revised 1 see-saw from sigma_A)", n));
            sqrt_start_rev_1 = tic;
            if sum == 1
                [sigma_A_rev_1_n, sigma_B_rev_1_n, beta_sqrt_rev_1(n), ...
                    sigma_A_sqrt_rev_1_n, sigma_B_sqrt_rev_1_n, ...
                    gamma_sqrt_rev_1(n), ...
                    gamma_sqrt_stop_rev_1(n)] = computesqrt_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    "sigma_A", "rev_1");
            else
                [sigma_A_rev_1_n, sigma_B_rev_1_n, beta_sqrt_rev_1(n), ...
                    sigma_A_sqrt_rev_1_n, sigma_B_sqrt_rev_1_n, ...
                    gamma_sqrt_rev_1(n), ...
                    gamma_sqrt_stop_rev_1(n)] = computesqrt(Pi, dim_A, ...
                                                    dim_B, threshold, ...
                                                    N_seeSaw, "sigma_A", ...
                                                    "rev_1");
            end
            if should_save == 1
                sigma_A_rev_1{n} = sigma_A_rev_1_n;
                sigma_B_rev_1{n} = sigma_B_rev_1_n;
                sigma_A_sqrt_rev_1{n} = sigma_A_sqrt_rev_1_n;
                sigma_B_sqrt_rev_1{n} = sigma_B_sqrt_rev_1_n;
            end
            sqrt_time_rev_1(n) = toc(sqrt_start_rev_1);
        end
        
        % revised 1 see-saw from sigma_B
        if ~any(ismember('sqrt_rev_1_B', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Square root (revised 1 see-saw from sigma_B)", n));
            sqrt_start_rev_1_B = tic;
            if sum == 1
                [sigma_A_rev_1_B_n, sigma_B_rev_1_B_n, ...
                    beta_sqrt_rev_1_B(n), sigma_A_sqrt_rev_1_B_n, ...
                    sigma_B_sqrt_rev_1_B_n, gamma_sqrt_rev_1_B(n), ...
                    gamma_sqrt_stop_rev_1_B(n)] = computesqrt_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    "sigma_B", "rev_1");
            else
                [sigma_A_rev_1_B_n, sigma_B_rev_1_B_n, ...
                    beta_sqrt_rev_1_B(n), sigma_A_sqrt_rev_1_B_n, ...
                    sigma_B_sqrt_rev_1_B_n, gamma_sqrt_rev_1_B(n), ...
                    gamma_sqrt_stop_rev_1_B(n)] = computesqrt(Pi, dim_A, ...
                                                    dim_B, threshold, ...
                                                    N_seeSaw, "sigma_B", ...
                                                    "rev_1");
            end
            if should_save == 1
                sigma_A_rev_1_B{n} = sigma_A_rev_1_B_n;
                sigma_B_rev_1_B{n} = sigma_B_rev_1_B_n;
                sigma_A_sqrt_rev_1_B{n} = sigma_A_sqrt_rev_1_B_n;
                sigma_B_sqrt_rev_1_B{n} = sigma_B_sqrt_rev_1_B_n;
            end
            sqrt_time_rev_1_B(n) = toc(sqrt_start_rev_1_B);
        end
        
        % revised 2 see-saw from sigma_A
        if ~any(ismember('sqrt_rev_2_A', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Square root (revised 2 see-saw from sigma_A)", n));
            sqrt_start_rev_2 = tic;
            if sum == 1
                [sigma_A_rev_2_n, sigma_B_rev_2_n, beta_sqrt_rev_2(n), ...
                    sigma_A_sqrt_rev_2_n, sigma_B_sqrt_rev_2_n, ...
                    gamma_sqrt_rev_2(n), ...
                    gamma_sqrt_stop_rev_2(n)] = computesqrt_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    "sigma_A", "rev_2");
            else
                [sigma_A_rev_2_n, sigma_B_rev_2_n, beta_sqrt_rev_2(n), ...
                    sigma_A_sqrt_rev_2_n, sigma_B_sqrt_rev_2_n, ...
                    gamma_sqrt_rev_2(n), ...
                    gamma_sqrt_stop_rev_2(n)] = computesqrt(Pi, dim_A, ...
                                                    dim_B, threshold, ...
                                                    N_seeSaw, "sigma_A", ...
                                                    "rev_2");
            end
            if should_save == 1
                sigma_A_rev_2{n} = sigma_A_rev_2_n;
                sigma_B_rev_2{n} = sigma_B_rev_2_n;
                sigma_A_sqrt_rev_2{n} = sigma_A_sqrt_rev_2_n;
                sigma_B_sqrt_rev_2{n} = sigma_B_sqrt_rev_2_n;
            end
            sqrt_time_rev_2(n) = toc(sqrt_start_rev_2);
        end
        
        % revised 2 see-saw from sigma_B
        if ~any(ismember('sqrt_rev_2_B', dont_do))
            waitbar(n/N_ops, h, ...
  sprintf("Instance %d: Square root (revised 2 see-saw from sigma_B)", n));
            sqrt_start_rev_2_B = tic;
            if sum == 1
                [sigma_A_rev_2_B_n, sigma_B_rev_2_B_n, ...
                    beta_sqrt_rev_2_B(n), sigma_A_sqrt_rev_2_B_n, ...
                    sigma_B_sqrt_rev_2_B_n, gamma_sqrt_rev_2_B(n), ...
                    gamma_sqrt_stop_rev_2_B(n)] = computesqrt_sum(K_list, ...
                                                    L_list, num_summands, ...
                                                    dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    "sigma_B", "rev_2");
            else
                [sigma_A_rev_2_B_n, sigma_B_rev_2_B_n, ...
                    beta_sqrt_rev_2_B(n), sigma_A_sqrt_rev_2_B_n, ...
                    sigma_B_sqrt_rev_2_B_n, gamma_sqrt_rev_2_B(n), ...
                    gamma_sqrt_stop_rev_2_B(n)] = computesqrt(Pi, dim_A, ...
                                                    dim_B, threshold, ...
                                                    N_seeSaw, "sigma_B", ...
                                                    "rev_2");
            end
            if should_save == 1
                sigma_A_rev_2_B{n} = sigma_A_rev_2_B_n;
                sigma_B_rev_2_B{n} = sigma_B_rev_2_B_n;
                sigma_A_sqrt_rev_2_B{n} = sigma_A_sqrt_rev_2_B_n;
                sigma_B_sqrt_rev_2_B{n} = sigma_B_sqrt_rev_2_B_n;
            end
            sqrt_time_rev_2_B(n) = toc(sqrt_start_rev_2_B);
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
        end
        
        %% gamma_rand: Random see-saw
        % simple see-saw
        if ~any(ismember('rand_simple', dont_do))
            waitbar(n/N_ops, h, ...
                    sprintf("Instance %d: Random (simple see-saw)", n));
            rand_start = tic;
            if sum == 1
                [sigma_A_rand_n, sigma_B_rand_n, gamma_rand(n), ...
                    gamma_rand_stop(n)] = computegamma_rand_sum(K_list, ...
                                            L_list, num_summands, dim_A, ...
                                            dim_B, threshold, N_seeSaw, ...
                                            N_rand, "simple");
            else
                [sigma_A_rand_n, sigma_B_rand_n, gamma_rand(n), ...
                    gamma_rand_stop(n)] = computegamma_rand(Pi, dim_A, ...
                                            dim_B, threshold, N_seeSaw, ...
                                                        N_rand, "simple");
            end
            if should_save == 1
                sigma_A_rand{n} = sigma_A_rand_n;
                sigma_B_rand{n} = sigma_B_rand_n;
            end
            rand_time(n) = toc(rand_start);
        end
        
        % revised 1 see-saw
        if ~any(ismember('rand_rev_1', dont_do))
            waitbar(n/N_ops, h, ...
                    sprintf("Instance %d: Random (revised 1 see-saw)", n));
            rand_start_rev_1 = tic;
            if sum == 1
                [sigma_A_rand_rev_1_n, sigma_B_rand_rev_1_n, ...
                    gamma_rand_rev_1(n), ...
                gamma_rand_stop_rev_1(n)] = computegamma_rand_sum(K_list, ...
                                                L_list, num_summands, ...
                                                dim_A, dim_B, threshold, ...
                                                N_seeSaw, N_rand, "rev_1");
            else
                [sigma_A_rand_rev_1_n, sigma_B_rand_rev_1_n, ...
                    gamma_rand_rev_1(n), ...
                    gamma_rand_stop_rev_1(n)] = computegamma_rand(Pi, ...
                                                    dim_A, dim_B, ...
                                                    threshold, N_seeSaw, ...
                                                    N_rand, "rev_1");
            end
            if should_save == 1
                sigma_A_rand_rev_1{n} = sigma_A_rand_rev_1_n;
                sigma_B_rand_rev_1{n} = sigma_B_rand_rev_1_n;
            end
            rand_time_rev_1(n) = toc(rand_start_rev_1);
        end
        
        % revised 2 see-saw
        if ~any(ismember('rand_rev_2', dont_do))
            waitbar(n/N_ops, h, ...
                    sprintf("Instance %d: Random (revised 2 see-saw)", n));
            rand_start_rev_2 = tic;
            if sum == 1
                [sigma_A_rand_rev_2_n, sigma_B_rand_rev_2_n, ...
                    gamma_rand_rev_2(n), ...
                gamma_rand_stop_rev_2(n)] = computegamma_rand_sum(K_list, ...
                                                L_list, num_summands, ...
                                                dim_A, dim_B, threshold, ...
                                                N_seeSaw, N_rand, "rev_2");
            else
                [sigma_A_rand_rev_2_n, sigma_B_rand_rev_2_n, ...
                    gamma_rand_rev_2(n), ...
                gamma_rand_stop_rev_2(n)] = computegamma_rand(Pi, dim_A, ...
                                                    dim_B, threshold, ...
                                                    N_seeSaw, N_rand, ...
                                                    "rev_2");
            end
            if should_save == 1
                sigma_A_rand_rev_2{n} = sigma_A_rand_rev_2_n;
                sigma_B_rand_rev_2{n} = sigma_B_rand_rev_2_n;
            end
            rand_time_rev_2(n) = toc(rand_start_rev_2);
        end
        
        %% beta: Smallest upper bound
        beta = nan;
        if ~any(ismember('PPT', dont_do)) || any(ismember('beta_PPT', ...
                var_list))
            beta = min(beta, beta_PPT(n));
        end
        if ~any(ismember('realignment', ...
                dont_do)) || any(ismember('beta_r', var_list))
            beta = min(beta, beta_r(n));
        end
        if ~any(ismember('SymExt', ...
                dont_do)) || any(ismember('beta_k', var_list))
            beta = min(beta, min(beta_k(:, n)));
        end
        if ~any(ismember('BosSymExt', ...
                dont_do)) || any(ismember('beta_prime_k', var_list))
            beta = min(beta, min(beta_prime_k(:, n)));
        end
        if strcmp(opt, "Pos")
            if ~any(ismember('sqrt_simple_A', ...
                    dont_do)) || any(ismember('beta_sqrt', var_list))
                beta = min(beta, beta_sqrt(n));
            end
            if ~any(ismember('sqrt_simple_B', ...
                    dont_do)) || any(ismember('beta_sqrt_B', var_list))
                beta = min(beta, beta_sqrt_B(n));
            end
            if ~any(ismember('sqrt_rev_1_A', ...
                    dont_do)) || any(ismember('beta_sqrt_rev_1', var_list))
                beta = min(beta, beta_sqrt_rev_1(n));
            end
            if ~any(ismember('sqrt_rev_1_B', ...
                    dont_do)) || any(ismember('beta_sqrt_rev_1_B', ...
                                                var_list))
                beta = min(beta, beta_sqrt_rev_1_B(n));
            end
            if ~any(ismember('sqrt_rev_2_A', ...
                    dont_do)) || any(ismember('beta_sqrt_rev_2', var_list))
                beta = min(beta, beta_sqrt_rev_2(n));
            end
            if ~any(ismember('sqrt_rev_2_B', ...
                    dont_do)) || any(ismember('beta_sqrt_rev_2_B', ...
                                                var_list))
                beta = min(beta, beta_sqrt_rev_2_B(n));
            end
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
        
        % Square root
        if ~any(ismember('sqrt_simple_A', ...
                dont_do)) || any(ismember('gamma_sqrt', var_list))
            gamma = max(gamma, gamma_sqrt(n));
        end
        if ~any(ismember('sqrt_simple_B', ...
                dont_do)) || any(ismember('gamma_sqrt_B', var_list))
            gamma = max(gamma, gamma_sqrt_B(n));
        end
        if ~any(ismember('sqrt_rev_1_A', ...
                dont_do)) || any(ismember('gamma_sqrt_rev_1', var_list))
            gamma = max(gamma, gamma_sqrt_rev_1(n));
        end
        if ~any(ismember('sqrt_rev_1_B', ...
                dont_do)) || any(ismember('gamma_sqrt_rev_1_B', var_list))
            gamma = max(gamma, gamma_sqrt_rev_1_B(n));
        end
        if ~any(ismember('sqrt_rev_2_A', ...
                dont_do)) || any(ismember('gamma_sqrt_rev_2', var_list))
            gamma = max(gamma, gamma_sqrt_rev_2(n));
        end
        if ~any(ismember('sqrt_rev_2_B', ...
                dont_do)) || any(ismember('gamma_sqrt_rev_2_B', var_list))
            gamma = max(gamma, gamma_sqrt_rev_2_B(n));
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
        end
    end
    if ~any(ismember('PPT', dont_do))
        if should_save == 1
            save(filename, 'rho_PPT', '-append', '-v7.3');
        end
        save(filename, 'beta_PPT', 'beta_PPT_time', '-append', '-v7.3');
    end
    if ~any(ismember('realignment', dont_do))
        if should_save == 1
            save(filename, 'rho_r', '-append', '-v7.3');
        end
        save(filename, 'beta_r', 'beta_r_time', '-append', '-v7.3');
    end
    if ~any(ismember('SymExt', dont_do))
        if should_save == 1
            save(filename, 'rho_k', '-append', '-v7.3');
        end
        save(filename, 'beta_k', 'beta_k_time', '-append', '-v7.3');
    end
    if ~any(ismember('BosSymExt', dont_do))
        if should_save == 1
            save(filename, 'rho_prime_k', '-append', '-v7.3');
        end
        save(filename, 'beta_prime_k', 'beta_prime_k_time', '-append', ...
                '-v7.3');
    end
    
    % Square root
    if ~any(ismember('sqrt_simple_A', dont_do))
        if should_save == 1
            save(filename, 'sigma_A', 'sigma_B', 'sigma_A_sqrt', ...
                    'sigma_B_sqrt', '-append', '-v7.3');
        end
        save(filename, 'beta_sqrt', 'gamma_sqrt', 'gamma_sqrt_stop', ...
                'sqrt_time', '-append', '-v7.3');
    end
    if ~any(ismember('sqrt_simple_B', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_start_B', 'sigma_B_start_B', ...
                 'sigma_A_sqrt_B', 'sigma_B_sqrt_B', '-append', '-v7.3');
        end
        save(filename, 'beta_sqrt_B', 'gamma_sqrt_B', ...
                'gamma_sqrt_stop_B', 'sqrt_time_B', '-append', '-v7.3');
    end
    if ~any(ismember('sqrt_rev_1_A', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rev_1', 'sigma_B_rev_1', ...
                    'sigma_A_sqrt_rev_1', 'sigma_B_sqrt_rev_1', ...
                    '-append', '-v7.3');
        end
        save(filename, 'beta_sqrt_rev_1', 'gamma_sqrt_rev_1', ...
                'gamma_sqrt_stop_rev_1', 'sqrt_time_rev_1', '-append', ...
                '-v7.3');
    end
    if ~any(ismember('sqrt_rev_1_B', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rev_1_B', 'sigma_B_rev_1_B', ...
                    'sigma_A_sqrt_rev_1_B', 'sigma_B_sqrt_rev_1_B', ...
                    '-append', '-v7.3');
        end
        save(filename, 'beta_sqrt_rev_1_B', 'gamma_sqrt_rev_1_B', ...
                'gamma_sqrt_stop_rev_1_B', 'sqrt_time_rev_1_B', ...
                '-append', '-v7.3');
    end
    if ~any(ismember('sqrt_rev_2_A', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rev_2', 'sigma_B_rev_2', ...
                    'sigma_A_sqrt_rev_2', 'sigma_B_sqrt_rev_2', ...
                    '-append', '-v7.3');
        end
        save(filename, 'beta_sqrt_rev_2', 'gamma_sqrt_rev_2', ...
                'gamma_sqrt_stop_rev_2', 'sqrt_time_rev_2', '-append', ...
                '-v7.3');
    end
    if ~any(ismember('sqrt_rev_2_B', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rev_2_B', 'sigma_B_rev_2_B', ...
                    'sigma_A_sqrt_rev_2_B', 'sigma_B_sqrt_rev_2_B', ...
                    '-append', '-v7.3');
        end
        save(filename, 'beta_sqrt_rev_2_B', 'gamma_sqrt_rev_2_B', ...
                'gamma_sqrt_stop_rev_2_B', 'sqrt_time_rev_2_B', ...
                '-append', '-v7.3');
    end
    
    % Maximally-mixed state
    if ~any(ismember('MM_simple', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_MM', 'sigma_B_MM', '-append', '-v7.3');
        end
        save(filename, 'gamma_MM', 'gamma_MM_stop', 'MM_time', ...
                '-append', '-v7.3');
    end
    if ~any(ismember('MM_rev_1', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_MM_rev_1', 'sigma_B_MM_rev_1', ...
                    '-append', '-v7.3');
        end
        save(filename, 'gamma_MM_rev_1', 'gamma_MM_stop_rev_1', ...
                'MM_time_rev_1', '-append', '-v7.3');
    end
    if ~any(ismember('MM_rev_2', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_MM_rev_2', 'sigma_B_MM_rev_2', ...
                    '-append', '-v7.3');
        end
        save(filename, 'gamma_MM_rev_2', 'gamma_MM_stop_rev_2', ...
                'MM_time_rev_2', '-append', '-v7.3');
    end
    
    % Uniform superposition state
    if ~any(ismember('US_simple', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_US', 'sigma_B_US', '-append', '-v7.3');
        end
        save(filename, 'gamma_US', 'gamma_US_stop', 'US_time', ...
                '-append', '-v7.3');
    end
    if ~any(ismember('US_rev_1', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_US_rev_1', 'sigma_B_US_rev_1', ...
                    '-append', '-v7.3');
        end
        save(filename, 'gamma_US_rev_1', 'gamma_US_stop_rev_1', ...
                'US_time_rev_1', '-append', '-v7.3');
    end
    if ~any(ismember('US_rev_2', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_US_rev_2', 'sigma_B_US_rev_2', ...
                    '-append', '-v7.3');
        end
        save(filename, 'gamma_US_rev_2', 'gamma_US_stop_rev_2', ...
                'US_time_rev_2', '-append', '-v7.3');
    end
    
    % Random starting point
    if ~any(ismember('rand_simple', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rand', 'sigma_B_rand', '-append', ...
                '-v7.3');
        end
        save(filename, 'gamma_rand', 'gamma_rand_stop', 'rand_time', ...
                '-append', '-v7.3');
    end
    if ~any(ismember('rand_rev_1', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rand_rev_1', 'sigma_B_rand_rev_1', ...
                    '-append', '-v7.3');
        end
        save(filename, 'gamma_rand_rev_1', 'gamma_rand_stop_rev_1', ...
                'rand_time_rev_1', '-append', '-v7.3');
    end
    if ~any(ismember('rand_rev_2', dont_do))
        if should_save == 1
            save(filename, 'sigma_A_rand_rev_2', 'sigma_B_rand_rev_2', ...
                    '-append', '-v7.3');
        end
        save(filename, 'gamma_rand_rev_2', 'gamma_rand_stop_rev_2', ...
                'rand_time_rev_2', '-append', '-v7.3');
    end
    
    save(filename, 'ToB', '-append', '-v7.3');
    
    close(h);
end
