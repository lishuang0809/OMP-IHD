function DCT_dict = create_overcomplete_dct_SL(N, L)
    % Generate an overcomplete DCT dictionary of size (N x L)
    % N: signal dimension (number of rows)
    % M: number of dictionary atoms (columns)

    % Coded by Shuang Li, Feb. 2025.


    % Create standard DCT basis
    DCT_basis = dctmtx(N); % Generates a square N x N DCT matrix

    % Expand the basis to form an overcomplete dictionary
    % Here we tile different frequency components by repeating the basis
    DCT_dict = zeros(N, L);

    for k = 1:L
        idx = mod(k-1, N) + 1; % Cycle through the DCT basis columns
        DCT_dict(:, k) = DCT_basis(:, idx); % Assign basis column
    end

    % Normalize columns to have unit norm
    DCT_dict = normc(DCT_dict);
end