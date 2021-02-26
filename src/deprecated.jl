
# deprecated formulation of the problem
module OldBasicMD

using LinearAlgebra

"""
    Boltzmann_likelihood(x, V, β)
Compute the unnormalized Boltzmann density, exp(-β V(x))
for potential V at inverse temperature β
"""
function Boltzmann_likelihood(x, V, β)
    w = exp(-β * V(x));
    return w
end

"""
    MALA_likelihood(X₀, X₁, gradV0, β, Δt)
Compute the likelihood of the MALA proposal X₀→X₁ for
potential with gradient gradV0 at inverse temperature β
and time step Δt
"""
function MALA_likelihood(X₀, X₁, gradV0, β, Δt)
    #norm of increment squared
    norminc2 = 0.0;
    for i = 1:length(X₀)
        norminc2 += (X₁[i]-X₀[i] + Δt * gradV0[i])^2;
    end
    return exp( -β * (norminc2) / (4*Δt) );
end

"""
    RWM(x₀, V, β, Δt, n_iters, return_trajectory=true)
Run an RWM sampler for distribution exp(-β V).  Proposal variance is
sqrt(2*Δt/β).  If `return_trajectory=true`, then the entire time series is
returned.  The acceptance rates are also returned.
"""
function RWM(x₀, V, β, Δt, n_iters; return_trajectory=true)

    # preallocate data structures
    X₀ = copy(x₀);
    Xp = similar(x₀);
    d = length(x₀)

    naccept = 0;

    if(return_trajectory)
        Xvals =zeros(d,n_iters);
        avals =zeros(n_iters);
    end

    V0 = V(X₀);
    gaussian_coef =  sqrt(2 * Δt /β);
    for j = 1:n_iters

        @. Xp = X₀ + gaussian_coef * randn();

        Vp = V(Xp);
        a = min(1, exp(β*(V0-Vp)))

        if rand()<a
            naccept = naccept+1;
            @. X₀ = Xp;
            V0 = Vp;
        end

        if(return_trajectory)
            @. Xvals[:,j] = X₀;
            avals[j] = naccept/j;
        end

    end
    if return_trajectory
        return Xvals, avals
    else
        return X₀, naccept/n_iters
    end

end # end RWM

"""
    RWM!(X₀, V, β, Δt, n_iters)
Run an in place RWM sampler for distribution exp(-β V).  Proposal variance is
sqrt(2*Δt/β).
"""
function RWM!(X₀, V, β, Δt, n_iters)

    # preallocate data structures
    Xp = similar(X₀);

    V0 = V(X₀);
    gaussian_coef =  sqrt(2 * Δt/β);
    for j = 1:n_iters

        @. Xp = X₀ + gaussian_coef * randn();

        Vp = V(Xp);
        a = min(1, exp(β*(V0-Vp)))

        if rand()<a
            @. X₀ = Xp;
            V0 = Vp;
        end
    end
    X₀
end # end RWM!

"""
    MALA(x₀, V, gradV!, β, Δt, n_iters, return_trajectory=true)
Run a MALA sampler for distribution exp(-β V).  Proposal variance is
sqrt(2*Δt/β).  If `return_trajectory=true`, then the entire time series is
returned.  The acceptance rates are also returned.
"""
function MALA(x₀, V, gradV!, β, Δt, n_iters; return_trajectory=true)

    # preallocate data structures
    X₀ = copy(x₀);
    Xp = similar(x₀);
    gradV0 = similar(x₀);
    gradVp = similar(x₀);
    d = length(x₀)

    naccept = 0;

    if(return_trajectory)
        Xvals =zeros(d,n_iters);
        avals =zeros(n_iters);
    end

    V0 = V(X₀);
    gaussian_coef =  sqrt(2 * Δt/β);
    gradV!(gradV0,X₀);
    for j = 1:n_iters

        @. Xp = X₀ - Δt * gradV0 + gaussian_coef * randn();

        Vp = V(Xp);

        gradV!(gradVp,Xp);

        g0 = MALA_likelihood(X₀, Xp, gradV0, β, Δt);
        gp = MALA_likelihood(Xp, X₀, gradVp, β, Δt);

        a = min(1, gp/g0 * exp(β*(V0-Vp)));

        if rand()<a
            naccept = naccept+1;
            @. X₀ = Xp;
            @. gradV0 = gradVp;
            V0 = Vp;
        end
        if(return_trajectory)
            @. Xvals[:,j] = X₀;
            avals[j] = naccept/j;
        end

    end
    if return_trajectory
        return Xvals, avals
    else
        return X₀, naccept/n_iters
    end

end # end MALA

"""
    MALA!(x₀, V, gradV!, β, Δt, n_iters, return_trajectory=true)
Run an in place MALA sampler for distribution exp(-β V).  Proposal variance is
sqrt(2*Δt/β).
"""
function MALA!(X₀, V, gradV!, β, Δt, n_iters)

    # preallocate data structures
    Xp = similar(X₀);
    gradV0 = similar(X₀);
    gradVp = similar(X₀);

    V0 = V(X₀);
    gaussian_coef =  sqrt(2 * Δt/β);
    gradV!(gradV0,X₀);
    for j = 1:n_iters

        @. Xp = X₀ - Δt * gradV0 + gaussian_coef * randn();

        Vp = V(Xp);

        gradV!(gradVp,Xp);

        g0 = MALA_likelihood(X₀, Xp, gradV0, β, Δt);
        gp = MALA_likelihood(Xp, X₀, gradVp, β, Δt);

        a = min(1, gp/g0 * exp(β*(V0-Vp)));

        if rand()<a
            @. X₀ = Xp;
            @. gradV0 = gradVp;
            V0 = Vp;
        end
    end
    X₀
end # end MALA!

"""
    EM(x₀, V, gradV!, β, Δt, n_iters, return_trajectory=true)
Run the Euler-Maruyama integrator for dX = -∇V(X)dt + sqrt(2/β)dW.   If
If `return_trajectory=true`, then the entire time series is returned.
"""
function EM(x₀, gradV!, β, Δt, n_iters; return_trajectory=true)

    # preallocate data structures
    X₀ = copy(x₀);
    gradV0 = similar(x₀);
    d = length(x₀)

    if(return_trajectory)
        Xvals =zeros(d, n_iters);
    end

    gaussian_coef = sqrt(2 * Δt / β);

    for j = 1:n_iters
        gradV!(gradV0,X₀);

        @. X₀ = X₀ - gradV0 * Δt +  gaussian_coef * randn();

        if(return_trajectory)
            @. Xvals[:,j] = X₀;
        end
    end
    if(return_trajectory)
        return Xvals
    else
        return X₀
    end

end

"""
    EM!(X₀, V, gradV!, β, Δt, n_iters, return_trajectory=true)
Run an in place Euler-Maruyama integrator for dX = -∇V(X)dt + sqrt(2/β)dW.
"""
function EM!(X₀, gradV!, β, Δt, n_iters)
    # preallocate data structures
    gradV0 = similar(X₀);

    gaussian_coef = sqrt(2 * Δt / β);

    for j = 1:n_iters
        gradV!(gradV0,X₀);

        @. X₀ = X₀ - gradV0 * Δt +  gaussian_coef * randn();
    end
    X₀
end # end EM!

"""
    BBK(q₀,p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)
Run the BBK integrator for interial Langevin   If `return_trajectory=true`,
then the entire time series is returned.
"""
function BBK(q₀, p₀, gradV!, β, γ, M, Δt, n_iters; return_trajectory=true)

    # diffusion coefficient
    σ=sqrt(2*γ/β);
    d = length(q₀);

    c₀ = sqrt(0.5*Δt);

    Q₀ = copy(q₀);
    P₀ = copy(p₀);
    Phalf = similar(p₀);
    gradV0 = similar(q₀);

    if return_trajectory
        Qvals = zeros(d,n_iters);
        Pvals = zeros(d,n_iters);
    end

    gradV!(gradV0,Q₀);
    for j = 1:n_iters
        @. Phalf = P₀ - 0.5 * Δt * gradV0 - 0.5 * Δt * γ/M * P₀ + c₀ * σ * randn();
        @. Q₀ = Q₀ + Δt/M * Phalf;
        gradV!(gradV0, Q₀);
        @. P₀ = (Phalf -0.5 * Δt * gradV0 + c₀ * σ * randn())/(1 + 0.5 * Δt * γ/M);

        if return_trajectory
            @. Qvals[:,j] = Q₀;
            @. Pvals[:,j] = P₀;
        end

    end
    if return_trajectory
        return Qvals, Pvals
    else
        return Q₀, P₀
    end
end # end BBK

"""
    BBK!(q₀,p₀, gradV!, β, Δt, n_iters)
Run an in place BBK integrator for interial Langevin.
"""
function BBK!(Q₀, P₀,gradV!, β, γ, M, Δt, n_iters)

    # diffusion coefficient
    σ=sqrt(2*γ/β);
    d = length(Q₀);

    Phalf = similar(P₀);
    gradV0 = similar(Q₀);

    gradV!(gradV0,Q₀);
    for j = 1:n_iters
        @. Phalf = P₀ - 0.5 * Δt * gradV0 - 0.5 * Δt * γ/M * P₀ + sqrt(0.5*Δt) * σ * randn();
        @. Q₀ = Q₀ + Δt/M * Phalf;
        gradV!(gradV0, Q₀);
        @. P₀ = (Phalf -0.5 * Δt * gradV0 + sqrt(0.5 * Δt) * σ * randn())/(1 + 0.5 * Δt * γ/M);

    end
    Q₀, P₀

end # end BBK!

"""
    ABOBA(q₀,p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)
Run the ABOBA integrator for interial Langevin.  If `return_trajectory=true`,
then the entire time series is returned.
"""
function ABOBA(q₀, p₀, gradV!, β, γ, M, Δt, n_iters; return_trajectory=true)

    c₀ = exp(-Δt * γ);
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β);
    d = length(q₀);

    Q₀ = copy(q₀);
    P₀ = copy(p₀);
    Qhalf = similar(q₀);
    Phalf = similar(p₀);
    P̂half = similar(p₀);
    gradVhalf = similar(q₀);
    sqrtM = sqrt.(M);

    if return_trajectory
        Qvals = zeros(d,n_iters);
        Pvals = zeros(d,n_iters);
    end

    for j in 1:n_iters
        @. Qhalf = Q₀ + 0.5 * Δt * P₀/M;
        gradV!(gradVhalf, Qhalf);
        @. Phalf = P₀ - 0.5 * Δt * gradVhalf;
        @. P̂half = c₀ * Phalf + c₁ * sqrtM * randn();
        @. P₀ = P̂half - 0.5 * Δt * gradVhalf;
        @. Q₀ = Qhalf + 0.5 * Δt * P₀/M;

        if return_trajectory
            @. Qvals[:,j] = Q₀;
            @. Pvals[:,j] = P₀;
        end

    end
    if return_trajectory
        return Qvals, Pvals
    else
        return Q₀, P₀
    end
end

"""
    ABOBA!(Q₀,P₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)
Run an in place ABOBA integrator for interial Langevin.
"""
function ABOBA!(Q₀, P₀, gradV!, β, γ, M, Δt, n_iters)

    c₀ = exp(-Δt * γ);
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β);
    d = length(Q₀);

    Qhalf = similar(Q₀);
    Phalf = similar(P₀);
    P̂half = similar(P₀);
    gradVhalf = similar(Q₀);
    sqrtM = sqrt.(M);

    for i in 1:n_iters
        @. Qhalf = Q₀ + 0.5 * Δt * P₀/M;
        gradV!(gradVhalf, Qhalf);
        @. Phalf = P₀ - 0.5 * Δt * gradVhalf;
        @. P̂half = c₀ * Phalf + c₁ * sqrtM * randn();
        @. P₀ = P̂half - 0.5 * Δt * gradVhalf;
        @. Q₀ = Qhalf + 0.5 * Δt * P₀/M;

    end
    Q₀, P₀

end

"""
    BAOAB(q₀,p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)
Run the BAOAB integrator for interial Langevin.  If `return_trajectory=true`,
then the entire time series is returned.
"""
function BAOAB(q₀, p₀, gradV!, β, γ, M, Δt, n_iters; return_trajectory=true)

    c₀ = exp(-Δt * γ);
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β);
    d = length(q₀);

    Q₀ = copy(q₀);
    P₀ = copy(p₀);
    Qhalf = similar(q₀);
    Phalf = similar(p₀);
    P̂half = similar(p₀);
    gradV0 = similar(q₀);
    sqrtM = sqrt.(M);

    if return_trajectory
        Qvals = zeros(d,n_iters);
        Pvals = zeros(d,n_iters);
    end

    gradV!(gradV0, Q₀);
    for j in 1:n_iters
        @. Phalf = P₀ - 0.5 * Δt * gradV0;
        @. Qhalf = Q₀ + 0.5 * Δt * Phalf/M;
        @. P̂half = c₀ * Phalf + c₁ * sqrtM * randn();
        @. Q₀ = Qhalf + 0.5 * Δt * P̂half/M;
        gradV!(gradV0, Q₀);
        @. P₀ = P̂half - 0.5 * Δt * gradV0;

        if return_trajectory
            @. Qvals[:,j] = Q₀;
            @. Pvals[:,j] = P₀;
        end

    end
    if return_trajectory
        return Qvals, Pvals
    else
        return Q₀, P₀
    end
end

"""
    BAOAB!(Q₀,P₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)
Run an in place BAOAB integrator for interial Langevin.
"""
function BAOAB!(Q₀, P₀, gradV!, β, γ, M, Δt, n_iters)

    c₀ = exp(-Δt * γ);
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β);
    d = length(Q₀);

    Qhalf = similar(Q₀);
    Phalf = similar(P₀);
    P̂half = similar(P₀);
    gradV0 = similar(Q₀);
    sqrtM = sqrt.(M);

    gradV!(gradV0, Q₀);
    for j in 1:n_iters
        @. Phalf = P₀ - 0.5 * Δt * gradV0;
        @. Qhalf = Q₀ + 0.5 * Δt * Phalf/M;
        @. P̂half = c₀ * Phalf + c₁ * sqrtM * randn();
        @. Q₀ = Qhalf + 0.5 * Δt * P̂half/M;
        gradV!(gradV0, Q₀);
        @. P₀ = P̂half - 0.5 * Δt * gradV0;

    end
    Q₀, P₀

end

"""
    Verlet(q₀,p₀, gradV!, M, Δt, n_iters, return_trajectory=true)
Run the Verlet integrator.  If `return_trajectory=true`, then the entire time
series is returned.
"""
function Verlet(q₀, p₀, gradV!, M, Δt, n_iters; return_trajectory=true)

    d = length(q₀);

    Q₀ = copy(q₀);
    P₀ = copy(p₀);
    gradV0 = similar(Q₀);
    Phalf = similar(P₀);

    if return_trajectory
        Qvals = zeros(d,n_iters);
        Pvals = zeros(d,n_iters);
    end

    gradV!(gradV0, Q₀);
    for j in 1:n_iters
        @. Phalf = P₀ - 0.5 * Δt * gradV0;
        @. Q₀ = Q₀ + Δt * Phalf/M;
        gradV!(gradV0,Q₀);
        @. P₀ = Phalf - 0.5 * Δt * gradV0;

        if return_trajectory
            @. Qvals[:,j] = Q₀;
            @. Pvals[:,j] = P₀;
        end

    end

    if return_trajectory
        return Qvals, Pvals
    else
        return Q₀, P₀
    end


end

"""
    Verlet!(q₀,p₀, gradV!, M, Δt, n_iters, return_trajectory=true)
Run an in place Verlet integrator.
"""
function Verlet!(Q₀,P₀, gradV!, M, Δt, n_iters)

    gradV0 = similar(Q₀);
    Phalf = similar(P₀);

    gradV!(gradV0, Q₀);
    for j in 1:n_iters
        @. Phalf = P₀ - 0.5 * Δt * gradV0;
        @. Q₀ = Q₀ + Δt * Phalf/M;
        gradV!(gradV0,Q₀);
        @. P₀ = Phalf - 0.5 * Δt * gradV0;

    end
    Q₀, P₀
end


"""
    HMC(x₀, V, gradV!, β, M, Δt, nΔt, n_iters)
Run the HMC sampler
"""

function HMC(x₀, V, gradV!, β, M, Δt, nΔt, n_iters; return_trajectory = true)

    Q₀ = copy(x₀);
    P₀ = similar(x₀);
    Qp = similar(x₀);
    Pp = similar(x₀);
    gradVp = similar(x₀);
    Phalf = similar(P₀);
    d = length(x₀);

    naccept = 0;

    if(return_trajectory)
        Xvals =zeros(d,n_iters);
        avals =zeros(n_iters);
    end

    V₀ = V(Q₀);
    for j in 1:n_iters
        @. P₀ = sqrt(M/β) * randn();
        E₀ = V₀ + 0.5 * P₀⋅ (P₀ ./M);

        # run Verlet integrator
        Qp = copy(Q₀);
        Pp = copy(P₀);
        gradV!(gradVp, Qp);
        for k in 1:nΔt
            @. Phalf = Pp - 0.5 * Δt * gradVp;
            @. Qp = Qp + Δt * Phalf/M;
            gradV!(gradVp,Qp);
            @. Pp = Phalf - 0.5 * Δt * gradVp;
        end

        # compute energy
        Vp = V(Qp);
        Ep = Vp + 0.5 * Pp⋅ (Pp ./M);

        # accept/reject
        a = min(1, exp( β *(E₀-Ep)));
        if rand()<a
            naccept = naccept+1;
            @. Q₀ = Qp;
            V₀ = Vp;
        end

        if(return_trajectory)
            @. Xvals[:,j] = Q₀;
            avals[j] = naccept/j;
        end
    end

    if return_trajectory
        return Xvals, avals
    else
        return Q₀, naccept/n_iters
    end
end


"""
    HMC!(Q₀, V, gradV!, β, M, Δt, nΔt, n_iters)
In place HMC sampler
"""

function HMC!(Q₀, V, gradV!, β, M, Δt, nΔt, n_iters)

    P₀ = similar(x₀);
    Qp = similar(x₀);
    Pp = similar(x₀);
    gradVp = similar(x₀);
    Phalf = similar(P₀);
    d = length(x₀);

    V₀ = V(Q₀);
    for j in 1:n_iters
        @. P₀ = sqrt(M/β) * randn();
        E₀ = V₀ + 0.5 * P₀⋅ (P₀ ./M);

        # run Verlet integrator
        Qp = copy(Q₀);
        Pp = copy(P₀);
        gradV!(gradVp, Qp);
        for k in 1:nΔt
            @. Phalf = Pp - 0.5 * Δt * gradVp;
            @. Qp = Qp + Δt * Phalf/M;
            gradV!(gradVp,Qp);
            @. Pp = Phalf - 0.5 * Δt * gradVp;
        end

        # compute energy
        Vp = V(Qp);
        Ep = Vp + 0.5 * Pp⋅ (Pp ./M);

        # accept/reject
        a = min(1, exp( β *(E₀-Ep)));
        if rand()<a
            naccept = naccept+1;
            @. Q₀ = Qp;
            V₀ = Vp;
        end
    end
    Q₀
end

end # end module
