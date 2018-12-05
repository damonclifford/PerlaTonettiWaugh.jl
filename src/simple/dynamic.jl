# Implementation of the simple model with time-varying objects.
# Dynamic calculations, defined for each time ∈ t.
function f(du,u,p,t)
    @unpack L_1, L_2, z, r, μ, g, υ, π_tilde, T, ξ, z = p
    # Validate upwind scheme direction.
    # μ + υ^2/2 - g(t) < 0 || @warn "Drift must be strictly negative at all times"
    r = r(t)
    g = g(t)
    # Carry out calculations.
    L = (r - g - ξ*(μ - g) - ξ^2 * υ^2/2)*I - (μ + ξ*υ^2 - g)*L_1 - υ^2/2 * L_2 # (eq:A.9)
    mul!(du, L, u)
    du .-= π_tilde.(t, z) # discretized system of ODE for v (eq:12)
end

# Constructor
function simpleODE(params, settings)
    # Unpack necessary objects.
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, T, g = settings
    M = length(z)
    # Discretize the operator.
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, ξ) # L_1_minus ≡ L_1 is the only one we use.
    # Calculate the stationary solution.
    r_T = r(T)
    g_T = g(T)
    π_tilde_T = π_tilde.(Ref(T), z)
    L_T = (r_T - g_T - ξ*(μ - g_T) - ξ^2 * υ^2/2)*I - (μ + ξ*υ^2 - g_T)*L_1_minus - υ^2/2 * L_2 # (eq:A.9)
    # Solution to the rescaled differential equation.
    v_T = L_T \ π_tilde_T
    # Bundle as before.
    p = (L_1 = L_1_minus, L_2 = L_2, z = z, g = g, r = r, υ = υ, π_tilde = π_tilde, T = T, μ = μ, ξ = ξ)
    return ODEProblem(f, v_T, (T, 0.0), p)
end

# Implementation of the simple model with time-varying objects, represented by DAE
# Dynamic calculations, defined for each time ∈ t.
function f!(resid,du,u,p,t)
    @unpack L_1, L_2, z, r, μ, g, υ, π_tilde, T, M, ξ, x, ω = p
    # Carry out calculations.
    v_t = u[1:M]
    g_t = u[M+1]
    L = (r(t) - g_t - ξ*(μ - g_t) - ξ^2 * υ^2/2)*I - (μ - g_t + υ^2*ξ)*L_1 - υ^2/2 * L_2 # (eq:A.9)
    resid[1:M] .= L * v_t - π_tilde.(t, z)
    resid[1:M] .-= du[1:M] # discretized system of ODE for v (eq:12)
    resid[M+1] = v_t[1] + x(t) - dot(ω, v_t) # value matching condition (eq:13)
end

# Constructor
function simpleDAE(params, settings)
    # Unpack necessary objects.
    @unpack μ, υ, θ, r, x, ξ, π_tilde = params
    @unpack z, T, g = settings
    M = length(z)
    # Quadrature weighting
    ω = ω_weights(z, θ, ξ)
    # Discretize the operator.
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, ξ) # L_1_minus ≡ L_1 is the only one we use.
    # Calculate the stationary solution.
    r_T = r(T)
    g_T = g(T)
    L_T = (r_T - g_T - ξ*(μ-g_T) - ξ^2 * υ^2/2)*I - (μ + ξ*υ^2 - g_T)*L_1_minus - υ^2/2 * L_2 # (eq:A.9)
    v_T = L_T \ π_tilde.(Ref(T), z)
    # Bundle as before.
    p = (L_1 = L_1_minus, L_2 = L_2, z = z, g = g, r = r, υ = υ, π_tilde = π_tilde, T = T, μ = μ, g_T = g_T, M = M, ξ = ξ, x = x, ω = ω)
    # Other objects
    u_T = [v_T; g_T]
    resid_M1 = zeros(M+1)
    return DAEProblem(f!, resid_M1, u_T, (T, 0.0), differential_vars = [fill(true, M); false], p)
end
