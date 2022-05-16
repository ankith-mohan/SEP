**SEP** is a MATLAB toolbox for computing bounds on the separability problem.

## Installation
1. Download [SEP](https://github.com/ankith-mohan/SEP/archive/refs/heads/main.zip)
2. Unzip this file in your MATLAB scripts directory
3. Download [QETLAB](http://www.qetlab.com/Main_Page)

## Getting started: Two introdutory examples
The first example will cover the generation of a 4 x 4 positive semidefinite matrix and computing the upper and lower bounds on this matrix.
    >> P = genMat(4, "Pos", 0)

    P =



    >> % Upper bounds
    >> [rho_PPT, beta_PPT] = computebeta_PPT(P)

    rho_PPT =



    beta_PPT =



    >> [rho_k, beta_k] = computebeta_k(P, 2, 2, 1)

    rho_k =



    beta_k =



    >> [rho_kr, beta_kr] = computebeta_kr(P, 2, 2, 1)

    rho_kr =



    beta_kr =



    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k(P, 2, 2, 1)

    rho_prime_k =



    beta_prime_k =



    >> [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr(P, 2, 2, 1)

    rho_prime_kr =



    beta_prime_kr =



    >> [rho_r, beta_r] = computebeta_r(P)

    rho_r =



    beta_r =



    >> % Lower bounds
    >> [sigma_A, sigma_B, gamma, ~] = computegamma_part(P, ones(2)/2, 2, 1e-13, 100, "sigma_A")

    sigma_A =



    sigma_B =



    gamma =



    >> [sigma_A_rev_1, sigma_B_rev_1, gamma_rev_1, ~] = computegamma_part_rev_1(P, ones(2)/2, 2, 1e-13, 100, "sigma_A")

    sigma_A_rev_1 =



    sigma_B_rev_1 =



    gamma_rev_1 =



    >> [sigma_A_rev_2, sigma_B_rev_2, gamma_rev_2, ~] = computegamma_part_rev_2(P, ones(2)/2, 2, 1e-13, 100, "sigma_A")

    sigma_A_rev_2 =



    sigma_B_rev_2 =



    gamma_rev_2 =



    >> [sigma_A_rand, sigma_B_rand, gamma_rand, ~] = computegamma_rand(P, 2, 2, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand =



    sigma_B_rand =



    gamma_rand =



The second example will describe the generation of $$\Pi = \sum\limits_{i=1}^{10} K_i \otimes L_i$$ where each $$K_i$$ and $$L_i$$ are 4 x 4 Hermitian matrices, and the computation of lower and upper bounds for this matrix.

    >> genOtimesHouseholder("Pi.mat", 4, 4, 10, "Herm")
    >> load("Pi.mat", "K_list", "L_list")
    >> % Upper bounds
    >> [rho_PPT, beta_PPT] = computebeta_PPT_sum(K_list, L_list, 10)
    
    rho_PPT =



    beta_PPT =



    >> [rho_k, beta_k] = computebeta_k_sum(K_list, L_list, 10, 4, 1, 1)

    rho_k =



    beta_k =



    >> [rho_kr, beta_kr] = computebeta_kr_sum(K_list, L_list, 10, 4, 1, 1)

    rho_kr =



    beta_kr =



    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k_sum(K_list, L_list, 10, 4, 1, 1)

    rho_prime_k =



    beta_prime_k =



    >> [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr_sum(K_list, L_list, 10, 4, 1, 1)

    rho_prime_kr =



    beta_prime_kr =



    >> [rho_r, beta_r] = computebeta_r_sum(K_list, L_list, 10)

    rho_r =



    beta_r =



    >> % Lower bounds
    >> [sigma_A, sigma_B, gamma, ~] = computegamma_part_sum(K_list, L_list, 10, ones(4)/4, 4, 1e-13, 100, "sigma_A")

    sigma_A =



    sigma_B =



    gamma =



    >> [sigma_A_rev_1, sigma_B_rev_1, gamma_rev_1, ~] = computegamma_part_rev_1_sum(K_list, L_list, 10, ones(4)/4, 4, 1e-13, 100, "sigma_A")

    sigma_A_rev_1 =



    sigma_B_rev_1 =



    gamma_rev_1 =



    >> [sigma_A_rev_2, sigma_B_rev_2, gamma_rev_2, ~] = computegamma_part_rev_2_sum(K_list, L_list, 10, ones(4)/4, 4, 1e-13, 100, "sigma_A")

    sigma_A_rev_2 =



    sigma_B_rev_2 =



    gamma_rev_2 =



    >> [sigma_A_rand, sigma_B_rand, gamma_rand, ~] = computegamma_rand_sum(K_list, L_list, 10, 4, 4, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand =



    sigma_B_rand =



    gamma_rand =



## [List of functions](funcs.md)