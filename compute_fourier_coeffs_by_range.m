function coeffs = compute_fourier_coeffs_by_range(vr, azimuths_rad, min_points)
% Computes the a0, a1, b1, a2, b2 Fourier coefficients at each range gate
% using the Browning & Wexler method.
%
% Inputs:
%   vr            - Radial velocity matrix [n_azimuths × n_ranges]
%   azimuths_rad  - Azimuth angles in radians [n_azimuths × 1]
%   min_points    - Minimum # of valid Vr values to attempt fit
%
% Output:
%   coeffs struct with fields a0, a1, b1, a2, b2 (1 × n_ranges)

    [~, n_rng] = size(vr);

    % Preallocate coefficient arrays
    a0 = NaN(1, n_rng);
    a1 = NaN(1, n_rng);
    b1 = NaN(1, n_rng);
    a2 = NaN(1, n_rng);
    b2 = NaN(1, n_rng);

    for r = 1:n_rng
        vr_col = vr(:, r);

        if all(isnan(vr_col))
            continue
        end

        valid = ~isnan(vr_col);
        if sum(valid) < min_points
            continue
        end

        beta = azimuths_rad(valid);   % [n_valid × 1]
        vr_valid = vr_col(valid);     % [n_valid × 1]

        % Design matrix: [1, cos(β), sin(β), cos(2β), sin(2β)]
        A = [ ...
            ones(size(beta)), ...
            cos(beta), ...
            sin(beta), ...
            cos(2 * beta), ...
            sin(2 * beta) ...
        ];

        % Solve least squares
        coeffs_array = A \ vr_valid;

        % Assign to arrays
        a0(r) = coeffs_array(1);
        a1(r) = coeffs_array(2);
        b1(r) = coeffs_array(3);
        a2(r) = coeffs_array(4);
        b2(r) = coeffs_array(5);
    end

    % Return as struct
    coeffs = struct( ...
        'a0', a0, ...
        'a1', a1, ...
        'b1', b1, ...
        'a2', a2, ...
        'b2', b2 ...
    );
end
