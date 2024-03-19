function [S_w, c, c_s1, c_s2] = solve_system_of_equations(S_w, S, c, c_s1, c_s2, k_rw, k_ro, p_c, u_w)
    % Solve the system of equations (8), (13), (19), (20) along with initial and boundary conditions
    
    % Define AD-OO variables for automatic differentiation
    S_w = initVariablesADI(S_w);
    c = initVariablesADI(c);
    c_s1 = initVariablesADI(c_s1);
    c_s2 = initVariablesADI(c_s2);
    
    % Define variables for equations
    t = 0; % Time variable (not used in this step)
    
    % Define boundary conditions
    bc_type = 'dir';
    bc_value = 1 - S_o_0;
    
    % Define pore volume
    Vp = sum(pore_volume(G, rock));

    % Define transport solver
    transSol = initTransportADI(G, 'c', c, 'bc', bc_type, 'bcvalue', bc_value);
    
    % Define system of equations
    eq1 = eqnWaterFractionADI(G, rock, S_w, 'phi', phi, 'bc', bc_type, 'bcvalue', bc_value);
    eq2 = eqnWaterTransportADI(G, rock, S_w, transSol, 'bc', bc_type, 'bcvalue', bc_value);
    eq3 = eqnNanoparticleRetentionADI(G, rock, c, u_w, gamma_d, gamma_e, u_c);
    eq4 = eqnPoreThroatBlockageADI(G, rock, c, u_w, gamma_pt);
    
    % Solve the system of equations
    [S_w, c, c_s1, c_s2] = solveEquationsADI({eq1, eq2, eq3, eq4}, t, 'bc', bc_type, 'bcvalue', bc_value);
    
    % Extract solution from AD-OO variables
    S_w = double(S_w);
    c = double(c);
    c_s1 = double(c_s1);
    c_s2 = double(c_s2);
end
