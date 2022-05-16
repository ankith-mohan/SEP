**SEP** is a MATLAB toolbox for computing bounds on the separability problem.

## Installation
1. Download [SEP](https://github.com/ankith-mohan/SEP/archive/refs/heads/main.zip)
2. Unzip this file in your MATLAB scripts directory
3. Download [QETLAB](http://www.qetlab.com/Main_Page)

## Getting started: Two introdutory examples
The first example will cover the generation of a 4 x 4 positive semidefinite matrix and computing the upper and lower bounds on this matrix.

    >> P = genMat(4, "Pos", 0)

    P =

        0.6127 + 0.0000i  -0.0499 + 0.0773i  -0.1744 - 0.1410i   0.1570 + 0.1253i
       -0.0499 - 0.0773i   0.7745 + 0.0000i   0.0436 - 0.0953i  -0.0039 + 0.1294i
       -0.1744 + 0.1410i   0.0436 + 0.0953i   0.5292 + 0.0000i   0.2143 - 0.1308i
        0.1570 - 0.1253i  -0.0039 - 0.1294i   0.2143 + 0.1308i   0.4182 + 0.0000i

    >> % Upper bounds
    >> [rho_PPT, beta_PPT] = computebeta_PPT(P)

    rho_PPT =

        0.4113 + 0.0000i  -0.0367 - 0.2006i  -0.3943 - 0.0741i  -0.0010 + 0.1989i
       -0.0367 + 0.2006i   0.1011 + 0.0000i   0.0713 - 0.1857i  -0.0969 - 0.0182i
       -0.3943 + 0.0741i   0.0713 + 0.1857i   0.3914 + 0.0000i  -0.0349 - 0.1909i
       -0.0010 - 0.1989i  -0.0969 + 0.0182i  -0.0349 + 0.1909i   0.0962 + 0.0000i

    beta_PPT =

        0.8309

    >> [rho_k, beta_k] = computebeta_k(P, 2, 2, 1)

    rho_k =

        0.4113 + 0.0000i  -0.0367 - 0.2006i  -0.3943 - 0.0741i  -0.0010 + 0.1989i
       -0.0367 + 0.2006i   0.1011 + 0.0000i   0.0713 - 0.1857i  -0.0969 - 0.0182i
       -0.3943 + 0.0741i   0.0713 + 0.1857i   0.3914 + 0.0000i  -0.0349 - 0.1909i
       -0.0010 - 0.1989i  -0.0969 + 0.0182i  -0.0349 + 0.1909i   0.0962 + 0.0000i

    beta_k =

        0.8309

    >> [rho_kr, beta_kr] = computebeta_kr(P, 2, 2, 1)

    rho_kr =

        0.4113 + 0.0000i  -0.0367 - 0.2006i  -0.3943 - 0.0741i  -0.0010 + 0.1989i
       -0.0367 + 0.2006i   0.1011 + 0.0000i   0.0713 - 0.1857i  -0.0969 - 0.0182i
       -0.3943 + 0.0741i   0.0713 + 0.1857i   0.3914 + 0.0000i  -0.0349 - 0.1909i
       -0.0010 - 0.1989i  -0.0969 + 0.0182i  -0.0349 + 0.1909i   0.0962 + 0.0000i

    beta_kr =

        0.8309

    >> [rho_prime_k, beta_prime_k] = computebeta_prime_k(P, 2, 2, 1)

    rho_prime_k =

        0.4113 + 0.0000i  -0.0367 - 0.2006i  -0.3943 - 0.0741i  -0.0010 + 0.1989i
       -0.0367 + 0.2006i   0.1011 + 0.0000i   0.0713 - 0.1857i  -0.0969 - 0.0182i
       -0.3943 + 0.0741i   0.0713 + 0.1857i   0.3914 + 0.0000i  -0.0349 - 0.1909i
       -0.0010 - 0.1989i  -0.0969 + 0.0182i  -0.0349 + 0.1909i   0.0962 + 0.0000i

    beta_prime_k =

        0.8309

    >> [rho_prime_kr, beta_prime_kr] = computebeta_prime_kr(P, 2, 2, 1)

    rho_prime_kr =

        0.4113 + 0.0000i  -0.0367 - 0.2006i  -0.3943 - 0.0741i  -0.0010 + 0.1989i
       -0.0367 + 0.2006i   0.1011 + 0.0000i   0.0713 - 0.1857i  -0.0969 - 0.0182i
       -0.3943 + 0.0741i   0.0713 + 0.1857i   0.3914 + 0.0000i  -0.0349 - 0.1909i
       -0.0010 - 0.1989i  -0.0969 + 0.0182i  -0.0349 + 0.1909i   0.0962 + 0.0000i

    beta_prime_kr =

        0.8309

    >> [rho_r, beta_r] = computebeta_r(P)

    rho_r =

        0.4113 + 0.0000i  -0.0367 - 0.2006i  -0.3943 - 0.0741i  -0.0010 + 0.1989i
       -0.0367 + 0.2006i   0.1011 + 0.0000i   0.0713 - 0.1857i  -0.0969 - 0.0182i
       -0.3943 + 0.0741i   0.0713 + 0.1857i   0.3914 + 0.0000i  -0.0349 - 0.1909i
       -0.0010 - 0.1989i  -0.0969 + 0.0182i  -0.0349 + 0.1909i   0.0962 + 0.0000i

    beta_r =

        0.8309

    >> % Lower bounds
    >> [sigma_A, sigma_B, gamma, ~] = computegamma_part(P, ones(2)/2, 2, 1e-13, 100, "sigma_A")

    sigma_A =

        0.9693 + 0.0000i   0.0057 + 0.1725i
        0.0057 - 0.1725i   0.0307 + 0.0000i

    sigma_B =

        0.0546 + 0.0000i  -0.1269 + 0.1885i
       -0.1269 - 0.1885i   0.9454 + 0.0000i

    gamma =

        0.8234

    >> [sigma_A_rev_1, sigma_B_rev_1, gamma_rev_1, ~] = computegamma_part_rev_1(P, ones(2)/2, 2, 1e-13, 100, "sigma_A")

    sigma_A_rev_1 =

        0.9693 + 0.0000i   0.0057 + 0.1725i
        0.0057 - 0.1725i   0.0307 + 0.0000i

    sigma_B_rev_1 =

        0.0546 + 0.0000i  -0.1269 + 0.1885i
       -0.1269 - 0.1885i   0.9454 + 0.0000i

    gamma_rev_1 =

        0.8234

    >> [sigma_A_rev_2, sigma_B_rev_2, gamma_rev_2, ~] = computegamma_part_rev_2(P, ones(2)/2, 2, 1e-13, 100, "sigma_A")

    sigma_A_rev_2 =

        0.1188 + 0.0000i  -0.2835 + 0.1560i
       -0.2835 - 0.1560i   0.8812 + 0.0000i

    sigma_B_rev_2 =

        0.5247 + 0.0000i   0.4650 - 0.1821i
        0.4650 + 0.1821i   0.4753 + 0.0000i

    gamma_rev_2 =

        0.7403

    >> [sigma_A_rand, sigma_B_rand, gamma_rand, ~] = computegamma_rand(P, 2, 2, 1e-13, 100, 100, "rev_2", 2)

    sigma_A_rand =

        0.9469 + 0.0000i  -0.2184 - 0.0505i
       -0.2184 + 0.0505i   0.0531 + 0.0000i

    sigma_B_rand =

        0.4645 + 0.0000i  -0.4626 + 0.1863i
       -0.4626 - 0.1863i   0.5355 + 0.0000i

    gamma_rand =

        0.8070

The second example will describe the generation of $$\Pi = \sum\limits_{i=1}^{10} K_i \otimes L_i$$ where each $$K_i$$ and $$L_i$$ are 4 x 4 Hermitian matrices, and the computation of lower and upper bounds for this matrix.

    >> genOtimesHouseholder("Pi.mat", 4, 4, 10, "Herm")
    >> load("Pi.mat", "K_list", "L_list")
    >> % Upper bounds
    >> [~, beta_PPT] = computebeta_PPT_sum(K_list, L_list, 10)
    
    beta_PPT =

        1.6567

    >> [~, beta_k] = computebeta_k_sum(K_list, L_list, 10, 4, 1, 1)

    beta_k =

        1.6567

    >> [~, beta_kr] = computebeta_kr_sum(K_list, L_list, 10, 4, 1, 1)

    beta_kr =

        1.6567

    >> [~, beta_prime_k] = computebeta_prime_k_sum(K_list, L_list, 10, 4, 1, 1)

    beta_prime_k =

        1.6567

    >> [~, beta_prime_kr] = computebeta_prime_kr_sum(K_list, L_list, 10, 4, 1, 1)

    beta_prime_kr =

        1.6567

    >> [~, beta_r] = computebeta_r_sum(K_list, L_list, 10)

    beta_r =

        1.6567

    >> % Lower bounds
    >> [sigma_A, sigma_B, gamma, ~] = computegamma_part_sum(K_list, L_list, 10, ones(4)/4, 4, 1e-13, 100, "sigma_A")

    sigma_A =

        0.6438 + 0.0000i   0.1070 - 0.0633i  -0.0707 + 0.1366i   0.3559 - 0.2521i
        0.1070 + 0.0633i   0.0240 + 0.0000i  -0.0252 + 0.0157i   0.0840 - 0.0069i
       -0.0707 - 0.1366i  -0.0252 - 0.0157i   0.0367 + 0.0000i  -0.0926 - 0.0478i
        0.3559 + 0.2521i   0.0840 + 0.0069i  -0.0926 + 0.0478i   0.2954 + 0.0000i

    sigma_B =

        0.1966 + 0.0000i  -0.2022 - 0.2503i   0.0510 + 0.1896i  -0.0092 + 0.1256i
       -0.2022 + 0.2503i   0.5266 + 0.0000i  -0.2939 - 0.1300i  -0.1505 - 0.1409i
        0.0510 - 0.1896i  -0.2939 + 0.1300i   0.1961 + 0.0000i   0.1188 + 0.0415i
       -0.0092 - 0.1256i  -0.1505 + 0.1409i   0.1188 - 0.0415i   0.0807 + 0.0000i

    gamma =

        1.6567

    >> [sigma_A_rev_1, sigma_B_rev_1, gamma_rev_1, ~] = computegamma_part_rev_1_sum(K_list, L_list, 10, ones(4)/4, 4, 1e-13, 100, "sigma_A")

    sigma_A_rev_1 =

        0.6438 + 0.0000i   0.1070 - 0.0633i  -0.0707 + 0.1366i   0.3559 - 0.2521i
        0.1070 + 0.0633i   0.0240 + 0.0000i  -0.0252 + 0.0157i   0.0840 - 0.0069i
       -0.0707 - 0.1366i  -0.0252 - 0.0157i   0.0367 + 0.0000i  -0.0926 - 0.0478i
        0.3559 + 0.2521i   0.0840 + 0.0069i  -0.0926 + 0.0478i   0.2954 + 0.0000i

    sigma_B_rev_1 =

        0.1966 + 0.0000i  -0.2022 - 0.2503i   0.0510 + 0.1896i  -0.0092 + 0.1256i
       -0.2022 + 0.2503i   0.5266 + 0.0000i  -0.2939 - 0.1300i  -0.1505 - 0.1409i
        0.0510 - 0.1896i  -0.2939 + 0.1300i   0.1961 + 0.0000i   0.1188 + 0.0415i
       -0.0092 - 0.1256i  -0.1505 + 0.1409i   0.1188 - 0.0415i   0.0807 + 0.0000i

    gamma_rev_1 =

        1.6567

    >> [sigma_A_rev_2, sigma_B_rev_2, gamma_rev_2, ~] = computegamma_part_rev_2_sum(K_list, L_list, 10, ones(4)/4, 4, 1e-13, 100, "sigma_A")

    sigma_A_rev_2 =

        0.6438 + 0.0000i   0.1070 - 0.0633i  -0.0707 + 0.1366i   0.3559 - 0.2521i
        0.1070 + 0.0633i   0.0240 + 0.0000i  -0.0252 + 0.0157i   0.0840 - 0.0069i
       -0.0707 - 0.1366i  -0.0252 - 0.0157i   0.0367 + 0.0000i  -0.0926 - 0.0478i
        0.3559 + 0.2521i   0.0840 + 0.0069i  -0.0926 + 0.0478i   0.2954 + 0.0000i

    sigma_B_rev_2 =

        0.1966 + 0.0000i  -0.2022 - 0.2503i   0.0510 + 0.1896i  -0.0092 + 0.1256i
       -0.2022 + 0.2503i   0.5266 + 0.0000i  -0.2939 - 0.1300i  -0.1505 - 0.1409i
        0.0510 - 0.1896i  -0.2939 + 0.1300i   0.1961 + 0.0000i   0.1188 + 0.0415i
       -0.0092 - 0.1256i  -0.1505 + 0.1409i   0.1188 - 0.0415i   0.0807 + 0.0000i

    gamma_rev_2 =

        1.6567

    >> [sigma_A_rand, sigma_B_rand, gamma_rand, ~] = computegamma_rand_sum(K_list, L_list, 10, 4, 4, 1e-13, 100, 100, "rev_1", 2)

    sigma_A_rand =

        0.6438 + 0.0000i   0.1070 - 0.0633i  -0.0707 + 0.1366i   0.3559 - 0.2521i
        0.1070 + 0.0633i   0.0240 + 0.0000i  -0.0252 + 0.0157i   0.0840 - 0.0069i
       -0.0707 - 0.1366i  -0.0252 - 0.0157i   0.0367 + 0.0000i  -0.0926 - 0.0478i
        0.3559 + 0.2521i   0.0840 + 0.0069i  -0.0926 + 0.0478i   0.2954 + 0.0000i

    sigma_B_rand =

        0.1966 + 0.0000i  -0.2022 - 0.2503i   0.0510 + 0.1896i  -0.0092 + 0.1256i
       -0.2022 + 0.2503i   0.5266 + 0.0000i  -0.2939 - 0.1300i  -0.1505 - 0.1409i
        0.0510 - 0.1896i  -0.2939 + 0.1300i   0.1961 + 0.0000i   0.1188 + 0.0415i
       -0.0092 - 0.1256i  -0.1505 + 0.1409i   0.1188 - 0.0415i   0.0807 + 0.0000i

    gamma_rand =

        1.6567

## [List of functions](funcs.md)