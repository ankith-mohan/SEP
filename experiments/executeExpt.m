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