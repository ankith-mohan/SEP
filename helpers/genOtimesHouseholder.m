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