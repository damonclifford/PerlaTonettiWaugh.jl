# User settings.
# State grids.
    grids = [   range(0.0, 5.0, length = 500), # none of these are front-loaded to the degree we use in production, which means numerical quantities will differ a bit from algebraic
                range(0.0, 7.0, length = 700),
                unique([range(0.0, 5.0, length = 500)' range(5.0, 7.0, length = 200)']), # irregular
            ]
# Overall parameters.
    params = @with_kw (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.01, d = 2.3701, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053, ζ_p = ζ, χ_p = χ)
    baseline = params()
# Solver settings.
    initial_values = [  [0.25, 3.0, 1.0],
                        [0.0190, 1.434969, 1.06517] # ~ equilibrium value for these parameters.
                     ]

# Execute tests.
@testset "Tests with Vanilla Parameters" begin
    # Compute algebraic solutions.
    algebraic_sols = [stationary_algebraic(baseline, settings_defaults(z_ex = grid, stationary_x0 = (x, y) -> iv)) for grid in grids, iv in initial_values]
    # Test algebraic equilibrium quantities.
    # Growth rate tests.
    algebraic_gs = (x -> x.g).(algebraic_sols)
    @test var(algebraic_gs) < 1e-6 # Tests that the solutions are similar to one another.
    @test algebraic_gs[1] ≈ 0.01900455065415125 # Tests proximity to true value.
    # Ω tests.
    algebraic_Ωs = (x -> x.Ω).(algebraic_sols)
    @test var(algebraic_Ωs) < 1e-6
    @test algebraic_Ωs[1] ≈ 1.0651775565541244
    # z_hat tests.
    algebraic_zs = (x -> x.z_hat).(algebraic_sols)
    @test var(algebraic_zs) < 1e-6
    @test all(algebraic_zs .≈ 1.434969541725385)

    # Compute numerical solutions.
    numerical_sols = [stationary_numerical(baseline, settings_defaults(z_ex = grid)) for grid in grids]
    # Test numerical equilibrium quantities.
    numerical_gs = (x -> x.g).(numerical_sols)
    @test var(numerical_gs) < 1e-6 # Tests that the solutions are similar to one another.
    # Ω tests.
    numerical_Ωs = (x -> x.Ω).(numerical_sols)
    @test var(numerical_Ωs) < 1e-6
    # z_hat tests.
    numerical_zs = (x -> x.z_hat).(numerical_sols)
    @test var(numerical_zs) < 1e-6
    # Numerical residuals tests.
    for i in 1:length(numerical_sols)
        # Get values to test.
        grid = grids[i]
        sol = numerical_sols[i]
        # Compute interim quantities.
        ω = ω_weights(grid, baseline.θ, baseline.σ-1)
        value_matching = sol.v_tilde[1] - dot(sol.v_tilde, ω) + sol.x
        free_entry = sol.v_tilde[1] - sol.x*(1-baseline.χ)/baseline.χ
        # Test.
        @test value_matching < 1e-8
        @test free_entry < 1e-8
        @test sol.L_tilde ≈ sol.L_tilde_a + sol.L_tilde_E + sol.L_tilde_x
    end
end

@testset "Tests with Perturbed Parameters" begin
    # Generate new parameters.
    newparams = params(σ = 4.0) # Arbitrary change.
    # Compute algebraic solutions.
    algebraic_sols = [stationary_algebraic(newparams, settings_defaults(z_ex = grid, stationary_x0 = (x, y) -> iv)) for grid in grids, iv in initial_values]
    # Test consistency.
    algebraic_gs = (x -> x.g).(algebraic_sols)
    algebraic_Ωs = (x -> x.Ω).(algebraic_sols)
    algebraic_zs = (x -> x.z_hat).(algebraic_sols)
    @test var(algebraic_gs) < 1e-6 # Tests that the solutions are similar to one another.
    @test var(algebraic_Ωs) < 1e-6
    @test var(algebraic_zs) < 1e-6

    # Compute numerical solutions.
    numerical_sols = [stationary_numerical(newparams, settings_defaults(z_ex = grid)) for grid in grids]
    # Test consistency.
    numerical_gs = (x -> x.g).(numerical_sols)
    numerical_Ωs = (x -> x.Ω).(numerical_sols)
    numerical_zs = (x -> x.z_hat).(numerical_sols)
    @test var(numerical_gs) < 1e-6 # Tests that the solutions are similar to one another.
    @test var(numerical_Ωs) < 1e-6
    @test var(numerical_zs) < 1e-6

    # Numerical residuals tests.
    for i in 1:length(numerical_sols)
        # Get values to test.
        grid = grids[i]
        sol = numerical_sols[i]
        # Compute interim quantities.
        ω = ω_weights(grid, baseline.θ, baseline.σ-1)
        value_matching = sol.v_tilde[1] - dot(sol.v_tilde, ω) + sol.x
        free_entry = sol.v_tilde[1] - sol.x*(1-baseline.χ)/baseline.χ
        # Test.
        @test value_matching < 1e-8
        @test free_entry < 1e-8
        @test sol.L_tilde ≈ sol.L_tilde_a + sol.L_tilde_E + sol.L_tilde_x
    end
end

@testset "Fixed g Method" begin
    # comparison to full method
    d_0 = 3.0426
    d_T =  2.83834
    params = parameter_defaults(γ = 1)
    settings = settings_defaults();
    params_0 = merge(params, (d = d_0,))
    params_T = merge(params, (d = d_T,))

    full_0 = stationary_algebraic(params_0, settings)
    full_T = stationary_algebraic(params_T, settings)
    g_0 = full_0.g
    g_T = full_T.g

    constrained_0 = stationary_algebraic_given_g(g_0, params_0, settings)
    @test constrained_0.z_hat ≈ full_0.z_hat
    @test constrained_0.Ω ≈ full_0.Ω
    constrained_T = stationary_algebraic_given_g(g_T, params_T, settings)
    @test constrained_T.z_hat ≈ full_T.z_hat
    @test constrained_T.Ω ≈ full_T.Ω

    # perturbation tests
    g_0 += 0.0001
    g_T -= 0.0001
    constrained_0_perturbed = stationary_algebraic_given_g(g_0, params_0, settings)
    constrained_T_perturbed = stationary_algebraic_given_g(g_T, params_T, settings)
    @test abs(constrained_0_perturbed.z_hat - constrained_0.z_hat) < 1e-3
    @test abs(constrained_0_perturbed.Ω - constrained_0.Ω) < 1e-3
    @test abs(constrained_T_perturbed.z_hat - constrained_T.z_hat) < 1e-3
    @test abs(constrained_T_perturbed.Ω - constrained_T.Ω) < 1e-3
end

@testset "Fixed g, Ω method" begin 
    d_0 = 3.0426
    d_T =  2.83834
    params = parameter_defaults(γ = 1)
    settings = settings_defaults();
    params_0 = merge(params, (d = d_0,))
    params_T = merge(params, (d = d_T,))

    full_0 = stationary_algebraic(params_0, settings)
    full_T = stationary_algebraic(params_T, settings)
    g_0 = full_0.g
    Ω_0 = full_0.Ω
    g_T = full_T.g
    Ω_T = full_T.Ω

    constrained_0 = stationary_algebraic_given_g_Ω(g_0, Ω_0, params_0, settings)
    @test constrained_0.z_hat ≈ full_0.z_hat atol = 1e-10
    @test constrained_0.L_tilde ≈ full_0.L_tilde atol = 1e-10 # arbitrary sanity check quantity
    constrained_T = stationary_algebraic_given_g(g_T, params_T, settings)
    @test constrained_T.z_hat ≈ full_T.z_hat atol = 1e-10
    @test constrained_T.L_tilde ≈ full_T.L_tilde atol = 1e-10 # arbitrary sanity check quantity

    # fudging the d slightly
    d_0_fudged = d_0 - 0.0001
    params_0_fudged = merge(params_0, (d = d_0_fudged,))
    constrained_0_fudged = stationary_algebraic_given_g_Ω(g_0, Ω_0, params_0_fudged, settings)
    @test constrained_0_fudged.z_hat ≈ constrained_0.z_hat atol = 1e-4
    @test constrained_0_fudged.L_tilde ≈ constrained_0.L_tilde atol = 1e-4
end

@testset "fixed c method" begin 
    d_0 = 3.0426
    params = parameter_defaults(γ = 1)
    params_0 = merge(params, (d = d_0,))
    settings = settings_defaults()
    full_0 = stationary_algebraic(params_0, settings)
    
    sol_c = steady_state_from_c(full_0.c, full_0.z_hat, full_0.Ω, params_0, settings)
    for k in keys(sol_c)
        sol_c[k] isa Number && @test sol_c[k] ≈ full_0[k]
    end

    sol_c_perturbed = steady_state_from_c(full_0.c + 0.001, full_0.z_hat, full_0.Ω, params_0, settings)
    for k in keys(sol_c_perturbed)
        sol_c_perturbed[k] isa Number && k != :U_bar && @test sol_c_perturbed[k] ≈ full_0[k] atol = 0.1
    end
end

# @testset "fixed lambda method" begin 
#     d_0 = 3.0426
#     params = parameter_defaults(γ = 1)
#     params_0 = merge(params, (d = d_0,))
#     settings = settings_defaults()
#     full_0 = stationary_algebraic(params_0, settings)
    
#     sol_λ = steady_state_from_λ(full_0.λ_ii, full_0.g, full_0.Ω, params_0, settings)
#     for k in keys(sol_λ)
#         sol_λ[k] isa Number && @test sol_λ[k] ≈ full_0[k]
#     end

#     sol_λ_perturbed = steady_state_from_λ(full_0.λ_ii + 0.001, full_0.g, full_0.Ω, params_0, settings)
#     for k in keys(sol_λ_perturbed)
#         sol_λ_perturbed[k] isa Number && k != :U_bar && @test sol_λ_perturbed[k] ≈ full_0[k] atol = 0.1
#     end
# end
