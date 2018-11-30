# Implements the algebraic stationary solution for the full model. Returns the equilibrium quantities (g, Ω, π) determined by equations H.15-H.17.


# Gives us the (full algebraic) stationary solution for a set of params and an initial x.
function stationary_algebraic(params, init_x = defaultiv(params); kwargs...)
    @assert params.υ > 0 && params.κ > 0 # Parameter validation
    sol = nlsolve(vals -> stationary_algebraic_aux(vals, params), init_x; inplace = false, kwargs...)
    converged(sol) || throw(sol)
    g, z_hat, Ω  = sol.zero
    @assert z_hat > 1 && Ω > 0 && g > 0 # Validate parameters.
    staticvalues = staticvals(sol.zero, params)
    return merge(staticvalues, (g = g, z_hat = z_hat, Ω = Ω))
end

# Gives us the residuals for a point x in state-space and a set of params.
function stationary_algebraic_aux(vals, params)
    # Grab values and intermediate quantities.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals(vals, params)
    g, z_hat, Ω = vals
    # Validate parameters.
    # Calculate and assign residuals.
        big_denom = ν*(θ + ν)*(θ - σ + 1) # Part of H.16
        denom_1 = a*(g - r) # Part of H.16
        num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # Part of H.16
        num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ -1) + (ν + σ - 1)*(θ + ν - σ + 1)) # Part of H.16
    return [x/π_min - a*(χ/(1-χ))*(σ + ν - 1)/ν, 1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν), π_min - (1- L_tilde)/((σ -1)*z_bar^(σ-1))]
end

function staticvals(vals, params)
    # Unpack x.
    g, z_hat, Ω = vals
    # Unpack params.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    # Compute interim quantities.
    F(z) = 1 - z^(-θ) # H.1
    r = ρ + γ*g + δ # H.6
    ν = (μ-g)/υ^2 + sqrt(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)) # H.3
    a = (r - g - (σ - 1)*(μ - g + (σ - 1)*υ^2/2))^(-1) # H.4
    b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν + σ - 1) # H.5
    S = θ * (g - μ - θ * υ^2 /2) # H.2
    L_tilde = Ω * ((N-1)*(1-F(z_hat))*κ + (1-η)*ζ*(S + δ/χ)) # H.7
    z_bar = (Ω * (θ/(1 + θ - σ) + (N-1)*(1-F(z_hat))*d^(1-σ)*(z_hat^(-1 + σ)*θ/(1 + θ - σ))))^((σ-1)^(-1)) # H.8
    w = σ^(-1)*z_bar # H.10
    x = ζ * (1- η + η * Theta / w) # H.11
    π_min = (d^(σ-1) * κ)/(z_hat^(σ-1)) # Inversion of H.9
    return (F = F, r = r, ν = ν, a = a, b = b, S = S, L_tilde = L_tilde, z_bar = z_bar, w = w, x = x, π_min = π_min)
end

# Default initial values
defaultiv(params) = [0.01, 2.0, 1.0]

# Numerical method.
function stationary_numerical(params, z, init_x = defaultiv(params); kwargs...)
    # Unpack params and settings.
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    @assert params.υ > 0 && params.κ > 0 # Parameter validation

    # Discretization objects and quadrature weights.
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operators.
    ω = ω_weights(z, θ, σ-1) # Get quadrature weights for the distribution on the rescaled grid.

    # Define the system of equations we're solving.
    function stationary_numerical_given_vals(vals)
        g, z_hat, Ω = vals
        @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals([g, z_hat, Ω], params) # Grab static values.
        r_tilde = r - g - 0 # g_w = 0 at steady state, equation B.31 (PDF)
        ρ_tilde = r_tilde - (σ - 1)*(μ - g + (σ-1)*(υ^2/2)) # B.21 (PDF)
        L_T = (ρ_tilde * I - (μ - g + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2) # Operator for the ξ-rescaled v function (23, PDF)
        i = z -> z >= log(z_hat) ? 1 : 0 # Indicator function for next equation.
        π_tilde = z -> π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (eq:32)
        v_tilde = L_T \ π_tilde.(z) # discretized system of ODE for v, where v'(T) = 0 (eq:24)

        #=
            System of equations to be solved.
        =#

        # Value-matching condition.
        value_matching = v_tilde[1] - dot(v_tilde, ω) + x # (eq:25)

        # Free-entry condition.
        free_entry = v_tilde[1] - x*(1-χ)/χ # (eq:27)

        # Adoption threshold.
        adoption_threshold = π_min - (1 - L_tilde)/((σ-1)*z_bar^(σ-1)) # (H.17) -- equivalent with (eq:26)

        return [value_matching, free_entry, adoption_threshold]
    end

    function f(x)
        resids = stationary_numerical_given_vals(x)
        return sum(resids .* resids)
    end

    function g!(G::Vector, x::Vector)
        ForwardDiff.gradient!(G, f, x)
    end

    function fg!(x::Vector, grad::Vector)
        if length(grad) > 0 # gradient of f(x)
            g!(grad, x)
        end
        f(x)
    end

    # define the optimization problem
    opt = Opt(:LD_LBFGS, 3) # 2 indicates the length of `x`
    lower_bounds!(opt, fill(0.0, 3)) # find `x` above 0
    min_objective!(opt, fg!) # specifies that optimization problem is on minimization
    xtol_rel!(opt, -Inf)
    xtol_abs!(opt, -Inf)
    ftol_rel!(opt, -Inf)
    ftol_abs!(opt, -Inf)

    # solve the optimization problem
    (minf,minx,ret) = NLopt.optimize(opt, [0.02; 18.94; 17.07])
    g_T, z_hat_T, Ω_T = minx
    # Grab static objects at steady-state.
    staticvalues = staticvals([g_T, z_hat_T, Ω_T], params) # Grab static values.
    @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvalues
    # Recreate the steady-state objects using the solution in g, z_hat, Ω.
    r_tilde = r - g_T - 0 # g_w = 0 at steady state (B.31, PDF)
    ρ_tilde = r_tilde - (σ - 1)*(μ - g_T + (σ-1)*(υ^2/2)) # (B.21, PDF)
    L_T = (ρ_tilde * I - (μ-g_T + (σ-1)*υ^2)*L_1_minus - υ^2/2 * L_2) # Operator for the ξ-rescaled v function (23, PDF)
    i = z -> z >= log(z_hat_T) ? 1 : 0 # Indicator function for next equation.
    π_tilde = z -> π_min * (1 + (N-1)*d^(1-σ)*i(z)) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # (eq:32)
    v_tilde = L_T \ π_tilde.(z) # discretized system of ODE for v, where v'(T) = 0 (eq:24)
    # Carry out welfare calculations.
    # λ_ii = 1/(1 + (N-1)*z_hat_T^(σ-1-θ)*d^(1-σ)) # (eq:H.21)
    # c_bar = (θ/(1-σ+θ))^(1/(σ-1))*(1-L_tilde)*Ω_T^(1/(σ-1))*λ_ii # (eq:B.54)
    # U_bar = ρ*log(c_bar) + g_T # (eq:43)

    return merge(staticvalues, (g = g_T, z_hat = z_hat_T, Ω = Ω_T, v_tilde = v_tilde))
end
