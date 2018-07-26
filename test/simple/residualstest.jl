# Global constants. 
    # State grid. 
    x_min = 0.0 # Rename to X for consistency. 
    x_max = 5.0
    M = 100
    x=linspace(x_min,x_max,M) # Since we only care about the grid. 

    # Time grid. 
    T = 10.0
    N = 10
    t = linspace(0.0, T, N)

    # Functional parameters. 
    π_func = (t, x) -> exp(x)
    ζ_func = t -> ζ_val # Not idiosyncratic, per equation (4)
    r_func = t -> r_val # Not idiosyncratic.

    # Constant parameters. 
    σ_val = 0.02
    α_val = 2.1
    ζ_val = 14.5
    r_val = 0.05
    γ_val = 0.005
    
    # Param generators and param NTs. 
    params = @with_kw (γ = γ_val, r = r_val, ζ = ζ_val, α = α_val, σ = σ_val) # Callable generator. 
    params_const = params()
    params_func = params(r = r_func, ζ = ζ_func) 

# Solutions. 
    # Solve for the numerical stationary g_T as test. 
    result_ns = stationary_numerical_simple(params_const, x)
    g_stationary = result_ns.g # This is the level. 

    # Test that this residual is close to 0. 
    ourDist = Truncated(Exponential(1/α_val), x[1], x[end]) 
    ω = irregulartrapezoidweights(x, ourDist)
    @test result_ns.v[1] + ζ_val - dot(ω, result_ns.v) ≈ 0 atol = 1e-10

    # Create the interpolation object of g
    g_vector = g_stationary + 0.01 * t
    g_int=LinInterp(t, g_vector)
    g_func = t -> g_int(t) # Not idiosyncratic. 

    # Calculate residuals. 
    resid= calculate_residuals(params_func, π_func, g_func, x, T)

    # Show the residual vector. 
    @show resid