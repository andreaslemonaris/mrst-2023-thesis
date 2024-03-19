function [phi, K, k_rw, k_ro, p_c, u_w] = update_values(S_w, c, c_s1, c_s2, phi, K, k_rw, k_ro, p_c, u_w)
    % Update values of phi, K, k_rw, k_ro, p_c, u_w based on the solution obtained
    
    % Update porosity
    delta_phi = c_s1 + c_s2;
    phi = phi_0 - delta_phi;
    
    % Update permeability
    f = 1 - gamma_f * c_s2;
    K = K_0 * ((1 - f) * k_f + f * phi / phi_0) .^ l;
    
    % Update relative permeabilities
    k_rw = k_rw_0 * S.^a;
    k_ro = k_ro_0 * (1 - S).^b;
    
    % Update capillary pressure
    p_c = b_w * (S + epsilon_1).^(-a_w) + b_o * (1 - S + epsilon_2).^(-a_o);
    
    % Update water velocity
    lambda_w = k_rw_0 * S_w.^a / mu_w;
    lambda_o = k_ro_0 * (1 - S_w).^b / mu_o;
    lambda_t = lambda_w + lambda_o;
    f_w = lambda_w ./ lambda_t;
    u_w = K .* lambda_w .* f_o .* ((dp_c / dz - delta_rho * g) / mu_w);
end
