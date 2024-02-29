var documenterSearchIndex = {"docs":
[{"location":"samplers/metropolis1/#Metropolis-Samplers","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"","category":"section"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"These samplers include a Metropolis-Hastings style step that ensure that the Boltzmann distribution mu(x) propto e^-beta V(x) is exactly targeted; there is no bias associated with, for instance, a finite time step Δt.","category":"page"},{"location":"samplers/metropolis1/#Zeroth-Order-Methods","page":"Metropolis Samplers","title":"Zeroth Order Methods","text":"","category":"section"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"These are samplers which do not require the gradienet of the potential, ∇V.","category":"page"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"    RWM(V::TV, β::TF, Δt::TF) where{TV, TF<:AbstractFloat}","category":"page"},{"location":"samplers/metropolis1/#BasicMD.RWM-Union{Tuple{TF}, Tuple{TV}, Tuple{TV, TF, TF}} where {TV, TF<:AbstractFloat}","page":"Metropolis Samplers","title":"BasicMD.RWM","text":"RWM(V, β, Δt)\n\nSet up the RWM sampler for Boltzmann.\n\nFields\n\nV     - Potential\nβ     - Inverse temperature\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"samplers/metropolis1/#First-Order-Methods","page":"Metropolis Samplers","title":"First Order Methods","text":"","category":"section"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"These are samplers which require the gradienet of the potential, ∇V, and are in the spirit of first order in time discretizations.","category":"page"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"    MALA(V::TV,∇V!::TGV, β::TF, Δt::TF) where{TV, TGV, TF<:AbstractFloat}","category":"page"},{"location":"samplers/metropolis1/#BasicMD.MALA-Union{Tuple{TF}, Tuple{TGV}, Tuple{TV}, Tuple{TV, TGV, TF, TF}} where {TV, TGV, TF<:AbstractFloat}","page":"Metropolis Samplers","title":"BasicMD.MALA","text":"MALA(V, ∇V!, β, Δt)\n\nSet up the MALA sampler for overdamped Langevin.\n\nFields\n\nV     - Potential\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"samplers/metropolis1/#Second-Order-Methods","page":"Metropolis Samplers","title":"Second Order Methods","text":"","category":"section"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"These are samplers which require the gradienet of the potential, ∇V, and are in the spirit of second order in time discretizations.","category":"page"},{"location":"samplers/metropolis1/","page":"Metropolis Samplers","title":"Metropolis Samplers","text":"    HMC(V::TV, ∇V!::TGV, β::TF, M::TM, Δt::TF, nΔt::Int) where {TV, TGV, TF<:AbstractFloat,TM}","category":"page"},{"location":"samplers/metropolis1/#BasicMD.HMC-Union{Tuple{TM}, Tuple{TF}, Tuple{TGV}, Tuple{TV}, Tuple{TV, TGV, TF, TM, TF, Int64}} where {TV, TGV, TF<:AbstractFloat, TM}","page":"Metropolis Samplers","title":"BasicMD.HMC","text":"HMC(V, ∇V!, β, M, Δt, nΔt)\n\nSet up the HMC sampler for Boltzmann.\n\nFields\n\nV     - Potential\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nM     - Mass matrix\nΔt    - Time step\nnΔt   - Number of time steps to use in each Verlet run\n\n\n\n\n\n","category":"method"},{"location":"samplers/nonmetropolis1/#Non-Metropolis-Samplers","page":"Non-Metropolis Samplers","title":"Non-Metropolis Samplers","text":"","category":"section"},{"location":"samplers/nonmetropolis1/","page":"Non-Metropolis Samplers","title":"Non-Metropolis Samplers","text":"These methods do not include a Metropolis-Hastings step, and, consequently, will sample from a distribution, mu_Delta t(x), which is a biased approximation of mu(x) propto e^-beta V(x).   This bias vanishes with Δt, and is often negligible in comparison to the statistical variance error.","category":"page"},{"location":"samplers/nonmetropolis1/#First-Order-Methods","page":"Non-Metropolis Samplers","title":"First Order Methods","text":"","category":"section"},{"location":"samplers/nonmetropolis1/","page":"Non-Metropolis Samplers","title":"Non-Metropolis Samplers","text":"These methods are in the spirit of first order in time discretizations.","category":"page"},{"location":"samplers/nonmetropolis1/","page":"Non-Metropolis Samplers","title":"Non-Metropolis Samplers","text":"    EM(∇V!::TGV, β::TF, Δt::TF) where{TGV, TF<:AbstractFloat}","category":"page"},{"location":"samplers/nonmetropolis1/#BasicMD.EM-Union{Tuple{TF}, Tuple{TGV}, Tuple{TGV, TF, TF}} where {TGV, TF<:AbstractFloat}","page":"Non-Metropolis Samplers","title":"BasicMD.EM","text":"EM(∇V!, β, γ, M, Δt)\n\nSet up the EM integrator for overdamped Langevin.\n\nFields\n\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"samplers/nonmetropolis1/#Second-Order-Methods","page":"Non-Metropolis Samplers","title":"Second Order Methods","text":"","category":"section"},{"location":"samplers/nonmetropolis1/","page":"Non-Metropolis Samplers","title":"Non-Metropolis Samplers","text":"These methods are in the spirit of second order in time discretizations.","category":"page"},{"location":"samplers/nonmetropolis1/","page":"Non-Metropolis Samplers","title":"Non-Metropolis Samplers","text":"    ABOBA(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}\n    BAOAB(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}\n    BBK(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}\n    GJF(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}","category":"page"},{"location":"samplers/nonmetropolis1/#BasicMD.ABOBA-Union{Tuple{TM}, Tuple{TF}, Tuple{TGV}, Tuple{TGV, TF, TF, TM, TF}} where {TGV, TF<:AbstractFloat, TM}","page":"Non-Metropolis Samplers","title":"BasicMD.ABOBA","text":"ABOBA(∇V!, β, γ, M, Δt)\n\nSet up the ABOBA integrator for inertial Langevin.\n\nFields\n\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nγ     - Damping Coefficient\nM     - Mass (either scalar or vector)\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"samplers/nonmetropolis1/#BasicMD.BAOAB-Union{Tuple{TM}, Tuple{TF}, Tuple{TGV}, Tuple{TGV, TF, TF, TM, TF}} where {TGV, TF<:AbstractFloat, TM}","page":"Non-Metropolis Samplers","title":"BasicMD.BAOAB","text":"BAOAB(∇V!, β, γ, M, Δt)\n\nSet up the BAOAB integrator for inertial Langevin.\n\nFields\n\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nγ     - Damping Coefficient\nM     - Mass (either scalar or vector)\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"samplers/nonmetropolis1/#BasicMD.BBK-Union{Tuple{TM}, Tuple{TF}, Tuple{TGV}, Tuple{TGV, TF, TF, TM, TF}} where {TGV, TF<:AbstractFloat, TM}","page":"Non-Metropolis Samplers","title":"BasicMD.BBK","text":"BBK(∇V!, β, γ, M, Δt)\n\nSet up the BBK integrator for inertial Langevin.\n\nFields\n\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nγ     - Damping Coefficient\nM     - Mass (either scalar or vector)\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"samplers/nonmetropolis1/#BasicMD.GJF-Union{Tuple{TM}, Tuple{TF}, Tuple{TGV}, Tuple{TGV, TF, TF, TM, TF}} where {TGV, TF<:AbstractFloat, TM}","page":"Non-Metropolis Samplers","title":"BasicMD.GJF","text":"GJF(∇V!, β, γ, M, Δt)\n\nSet up the G-JF integrator for inertial Langevin.\n\nFields\n\n∇V!   - In place gradient of the potential\nβ     - Inverse temperature\nγ     - Damping Coefficient\nM     - Mass (either scalar or vector)\nΔt    - Time step\n\n\n\n\n\n","category":"method"},{"location":"examples/sample_obs1/#Sampling-Observables","page":"Sampling Observables","title":"Sampling Observables","text":"","category":"section"},{"location":"examples/sample_obs1/","page":"Sampling Observables","title":"Sampling Observables","text":"Pages = [\"sample_obs1.md\"]","category":"page"},{"location":"examples/sample_obs1/","page":"Sampling Observables","title":"Sampling Observables","text":"These examples are performed using the 2D EntropicSwitch potential from TestLandscapes.jl","category":"page"},{"location":"examples/sample_obs1/#RWM-Example","page":"Sampling Observables","title":"RWM Example","text":"","category":"section"},{"location":"examples/sample_obs1/","page":"Sampling Observables","title":"Sampling Observables","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\nusing TestLandscapes\n\nV = x->EntropicSwitch(x);\n\nβ = 3.0;\nx0 = [-1.0, 0.0];\nseed = 100;\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = RWM(V, β, Δt);\n\n# construct observables\nf₁ = x-> x[1]^2; # second moment of coordinate 1\nf₂ = x-> x[2]^2; # second moment of coordinate 2\nf₃ = x-> V(x);   # energy\nobs = (f₁, f₂, f₃);\n\nRandom.seed!(100);\n\nobservable_samples = sample_observables(x0, sampler,obs, options=MDOptions(n_iters=n_iters));\n\nplot(1:n_iters, cumsum(observable_samples[1,:])./(1:n_iters), label=\"E[(X₁)²]\")\nplot!(1:n_iters, cumsum(observable_samples[2,:])./(1:n_iters), label=\"E[(X₂)²]\")\nplot!(1:n_iters, cumsum(observable_samples[3,:])./(1:n_iters), label=\"E[V(X)]\")\nxlabel!(\"Iterate\")\n","category":"page"},{"location":"examples/sample_obs1/#HMC-Example","page":"Sampling Observables","title":"HMC Example","text":"","category":"section"},{"location":"examples/sample_obs1/","page":"Sampling Observables","title":"Sampling Observables","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\nusing ForwardDiff\nusing TestLandscapes\n\nV = x->EntropicSwitch(x);\ngradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);\n\nβ = 3.0;\nx0 = [-1.0, 0.0];\nseed = 100;\nM = 1.;\nnΔt = 10^1; # number of Verlet steps per HMC iteration\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = HMC(V, gradV!, β, M, Δt, nΔt);\n\n# construct observables\nf₁ = x-> x[1]^2; # second moment of coordinate 1\nf₂ = x-> x[2]^2; # second moment of coordinate 2\nf₃ = x-> V(x);   # energy\nobs = (f₁, f₂, f₃);\n\nRandom.seed!(100);\n\nobservable_samples = sample_observables(x0, sampler,obs, options=MDOptions(n_iters=n_iters));\n\nplot(1:n_iters, cumsum(observable_samples[1,:])./(1:n_iters), label=\"E[(X₁)²]\")\nplot!(1:n_iters, cumsum(observable_samples[2,:])./(1:n_iters), label=\"E[(X₂)²]\")\nplot!(1:n_iters, cumsum(observable_samples[3,:])./(1:n_iters), label=\"E[V(X)]\")\nxlabel!(\"Iterate\")\n","category":"page"},{"location":"examples/sample_obs1/#ABOBA-Example","page":"Sampling Observables","title":"ABOBA Example","text":"","category":"section"},{"location":"examples/sample_obs1/","page":"Sampling Observables","title":"Sampling Observables","text":"note: Coordinates for Inertial Samplers\nFor inertial samplers, which approximate the underdamped Langevin equation, since the trajectory x(t) = (q(t) p(t)), it is essential what when you are evaluating observables, you extract the position or momenum, as appropriate.","category":"page"},{"location":"examples/sample_obs1/","page":"Sampling Observables","title":"Sampling Observables","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\nusing ForwardDiff\n\nusing TestLandscapes\n\nV = x->EntropicSwitch(x);\ngradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);\n\nβ = 3.0;\nγ = 1.;\nM = 1.;\nq0 = [-1.0, 0.0];\np0 = [0.0, 0.0]; # as ABOBA is an inertial sampler, we need to specify a momentum\nx0 = [copy(q0), copy(p0)];\nseed = 100;\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = ABOBA(gradV!, β, γ, M, Δt);\nopts = MDOptions(n_iters=n_iters);\n\n# construct observables\nf₁ = x-> x[1][1]^2; # second moment of coordinate 1\nf₂ = x-> x[1][2]^2; # second moment of coordinate 2\nf₃ = x-> V(x[1]);   # energy\nobs = (f₁, f₂, f₃);\n\nRandom.seed!(100);\nobservable_samples = sample_observables(x0, sampler,obs, options=MDOptions(n_iters=n_iters));\n\nplot(1:n_iters, cumsum(observable_samples[1,:])./(1:n_iters), label=\"E[(X₁)²]\")\nplot!(1:n_iters, cumsum(observable_samples[2,:])./(1:n_iters), label=\"E[(X₂)²]\")\nplot!(1:n_iters, cumsum(observable_samples[3,:])./(1:n_iters), label=\"E[V(X)]\")\nxlabel!(\"Iterate\")","category":"page"},{"location":"examples/sample_traj1/#Sampling-Trajectories","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"","category":"section"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"Pages = [\"sample_traj1.md\"]","category":"page"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"Examples of the sample_trajectory command with different samplers.","category":"page"},{"location":"examples/sample_traj1/#RWM-Example","page":"Sampling Trajectories","title":"RWM Example","text":"","category":"section"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\n\nfunction V(x)\n    return (x[1]^2 -1)^2\nend\n\nβ = 5.0;\nx0 = [-1.0];\nseed = 100;\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = RWM(V, β, Δt);\n\nRandom.seed!(100);\n\nxvals, avals = sample_trajectory(x0, sampler, options=MDOptions(n_iters=n_iters));\n\nhistogram([x_[1] for x_ in xvals],label=\"RWM Samples\",normalize=true)\nxlabel!(\"x\")\nylabel!(\"Frequency\")","category":"page"},{"location":"examples/sample_traj1/#HMC-Example","page":"Sampling Trajectories","title":"HMC Example","text":"","category":"section"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\nusing ForwardDiff\n\nfunction V(x)\n    return (x[1]^2 -1)^2\nend\n\ngradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);\n\n\nβ = 5.0;\nx0 = [-1.0];\nseed = 100;\nM = 1.;\nnΔt = 10^1; # number of Verlet steps per HMC iteration\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = HMC(V, gradV!, β, M, Δt, nΔt);\n\nRandom.seed!(100);\n\nxvals, avals = sample_trajectory(x0, sampler, options=MDOptions(n_iters=n_iters));\n\nhistogram([x_[1] for x_ in xvals],label=\"HMC Samples\",normalize=true)\nxlabel!(\"x\")\nylabel!(\"Frequency\")","category":"page"},{"location":"examples/sample_traj1/#Euler-Maruyama-Example","page":"Sampling Trajectories","title":"Euler-Maruyama Example","text":"","category":"section"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\nusing ForwardDiff\n\nfunction V(x)\n    return (x[1]^2 -1)^2\nend\n\ngradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);\n\n\nβ = 5.0;\nx0 = [-1.0];\nseed = 100;\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = EM(gradV!, β, Δt);\n\nRandom.seed!(100);\n\nxvals = sample_trajectory(x0, sampler, options=MDOptions(n_iters=n_iters));\n\nhistogram([x_[1] for x_ in xvals],label=\"EM Samples\",normalize=true)\nxlabel!(\"x\")\nylabel!(\"Frequency\")","category":"page"},{"location":"examples/sample_traj1/#ABOBA-Example","page":"Sampling Trajectories","title":"ABOBA Example","text":"","category":"section"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"note: Coordinates for Inertial Samplers\nFor inertial samplers, which approximate the underdamped Langevin equation, since the trajectory x(t) = (q(t) p(t)), it is neccessary to provide both an initial position and an initial momentum.","category":"page"},{"location":"examples/sample_traj1/","page":"Sampling Trajectories","title":"Sampling Trajectories","text":"using Plots\nusing Printf\nusing Random\nusing BasicMD\nusing ForwardDiff\n\nfunction V(x)\n    return (x[1]^2 -1)^2\nend\n\ngradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);\n\n\nβ = 5.0;\nγ = 1.;\nM = 1.;\nq0 = [-1.0];\np0 = [0.0]; # as ABOBA is an inertial sampler, we need to specify a momentum\nx0 = [copy(q0), copy(p0)];\nseed = 100;\nΔt = 1e-1;\nn_iters = 10^4; # number of samples\n\nsampler = ABOBA(gradV!, β, γ, M, Δt);\nopts = MDOptions(n_iters=n_iters);\n\nRandom.seed!(100);\nxvals = sample_trajectory(x0, sampler, options= opts);\n\nhistogram([x_[1][1] for x_ in xvals],label=\"ABOBA Samples\",normalize=true)\nxlabel!(\"x\")\nylabel!(\"Frequency\")","category":"page"},{"location":"sample1/#Creating-and-Using-Samplers","page":"Sampling","title":"Creating and Using Samplers","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Pages = [\"sample1.md\"]\nDepth = 4","category":"page"},{"location":"sample1/#Sampling-Trajectories-and-Boltzmann","page":"Sampling","title":"Sampling Trajectories and Boltzmann","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"At their core, all of the included methods sample, approximately, from Botlzmann distributinos of the type mu(x) propto e^-beta V(x), where the user must specify the potential, V(x), along with an appropriate inverse temperature, β.  ","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"The potential, V(x), must be formulated so that its input argument, x is an array, even if the problem is in mathbbR^1.  For the scalar double well potential, V(x) = (x^2-1)^2, this would be implemneted as:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"function V(x)\n    return (x[1]^2 -1)^2\nend","category":"page"},{"location":"sample1/#Constructing-the-Sampler","page":"Sampling","title":"Constructing the Sampler","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Having defined the potential, a sampler object must first be defined.  For instance,","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"rwm_sampler = RWM(V, β, Δt);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"constructs the random walk Metropolis (RWM) sampler for the Boltzmann disribution with `time step''Δt`.","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Other samplers require additional arguments.  For instance to use the HMC sampler, we would call","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"sampler = HMC(V, gradV!, β, M, Δt, nΔt);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"where the additional arguments are:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"gradV! -  in-place implementation of ∇V(x), \nM - mass matrix\nnΔt - number of Verlet steps of size Δt in each HMC iteration.\n","category":"page"},{"location":"sample1/#Sampling-a-Trajectory","page":"Sampling","title":"Sampling a Trajectory","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"    sample_trajectory(::Tx, ::S; options = MDOptions()) where {Tx,S<:BasicMD.MetropolisSampler}","category":"page"},{"location":"sample1/#BasicMD.sample_trajectory-Union{Tuple{S}, Tuple{Tx}, Tuple{Tx, S}} where {Tx, S<:BasicMD.MetropolisSampler}","page":"Sampling","title":"BasicMD.sample_trajectory","text":"sample_trajectory(x₀, sampler; options=MDOptions())\n\nRun the sampler starting at x₀.  Number of iterations and interval between saves are set using the options argument.  For Metropolis samplers, the running acceptance rates are also resturned.\n\nFields\n\nx₀         - Starting position for sampler\nsampler   - Desired sampler\n\nOptional Fields\n\noptions   - Sampling options, including number of iteration\nconstraints - Constraints on the trajectory\n\n\n\n\n\n","category":"method"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Having created the sampler strucutre and chosen an initial point, x0, we call sample_trajectory:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"xvals, avals = sampler_trajectory(x0, sampler);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"For a Metropolis sampler, like RWM, we return:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"xvals - the array of sampled points\navals - the running average of the acceptance rate For non-Metropolis","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"samplers, the avals argument is not returned. ","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Examples:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Pages = [\"examples/sample_traj1.md\"]","category":"page"},{"location":"sample1/#Controlling-Output","page":"Sampling","title":"Controlling Output","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"The sample_trajectory command will allocates and returns an array of samples.  The number of iterations and the sampling frequency is controlled through the options argument.  By default, if this is not specified, the sampler will be called for 10^4 iterations, and record every iteration.  To change this, we construct an MDOptions structure and pass that in.","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Additionally, one may be in the setting where we do not need to record the full trajectory, but merely the position at the terminal iterate.  This is handled with","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"    sample_trajectory!(::Tx, ::S; options = MDOptions()) where {Tx,S<:BasicMD.MetropolisSampler}","category":"page"},{"location":"sample1/#BasicMD.sample_trajectory!-Union{Tuple{S}, Tuple{Tx}, Tuple{Tx, S}} where {Tx, S<:BasicMD.MetropolisSampler}","page":"Sampling","title":"BasicMD.sample_trajectory!","text":"sample_trajectory!(x, sampler; options=MDOptions())\n\nIn place applciation of the sampler to x.  Number of iterations are set using the options argument.\n\nFields\n\nx         - Starting position for sampler, modified in place\nsampler   - Desired sampler\n\nOptional Fields\n\noptions   - Sampling options, including number of iteration\n\n\n\n\n\n","category":"method"},{"location":"sample1/#Sampling-Observables","page":"Sampling","title":"Sampling Observables","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Often, one is not interested in the full trajectory X_t, but rather the time series observables computed on the trajectory.  Given a function fmathbbR^dto mathbbR, it may be satisfactory to simply know f(X_t).  When d is high dimensional, only saving the observavble can cut computational cost and storage. Typically, this is done in order to estimate ergodic averages,","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"mathbbE_muf(X) = int f(x)mu(dx) approx frac1nsum_n=1^n f(X_t_n)","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"This is accomplished using","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"    sample_observables(x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}};\n    options = MDOptions()) where {Tx,S<:BasicMD.AbstractSampler,NO}","category":"page"},{"location":"sample1/#BasicMD.sample_observables-Union{Tuple{NO}, Tuple{S}, Tuple{Tx}, Tuple{Tx, S, Tuple{Vararg{Function, NO}}}} where {Tx, S<:BasicMD.AbstractSampler, NO}","page":"Sampling","title":"BasicMD.sample_observables","text":"sample_observables(x₀, sampler, observables; options=MDOptions())\n\nRun the sampler starting at x₀, evaluating the trajectory on a tuple of observables scalar functions.  Number of iterations and interval between saves are set using the options argument.  Only the computed observables are returned.\n\nFields\n\nx         - Starting position for sampler, modified in place\nsampler   - Desired sampler\nobservables - Observables on which to evaluate the trajectory\n\nOptional Fields\n\nTO- Observable data type, if needed, should be entered as the first argument\noptions   - Sampling options, including number of iteration\n\n\n\n\n\n","category":"method"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"This is quite similar to sample_trajectory, except that one must pass an additional argument, a structure containing the desired observables.  Additionally, only the time series of observables is returned, regardless of whether a Metropolis sampler is used or not.  ","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"The observables argument should be a tuple of scalar valued functions, i.e. for a single observable:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"f(x) = x[1]; # first component\nobs = (f,);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"and for multiple observables:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"f₁(x) = x[1]; # first component\nf₂(x) = V(x); # energy\nobs = (f₁,f₂);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"When calling the function,","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"observable_traj = sample_observables(x₀, sampler,obs, options=opts);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"the returned object, observable_traj is a matrix.  Each row corresponding to an individual observable, recorded at the times specified by the MDOptions.","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Examples:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Pages = [\"examples/sample_obs1.md\"]","category":"page"},{"location":"sample1/#Sampling-with-Constraints-(EXPERIMENTAL)","page":"Sampling","title":"Sampling with Constraints (EXPERIMENTAL)","text":"","category":"section"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"Elementary tools for enforcing constraints on the system, like X_t in A or g(X_t)=0 have been implemented.  This is accomplished by passing a constraints structure to one of sample_trajectory!, sample_trajectory, or sample_obsevables:","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"traj = sample_trajectory(x₀, sampler, constraints);","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"The constraints are constructed with","category":"page"},{"location":"sample1/","page":"Sampling","title":"Sampling","text":"    trivial_constraint!","category":"page"},{"location":"sample1/#BasicMD.trivial_constraint!","page":"Sampling","title":"BasicMD.trivial_constraint!","text":"trivial_constraint!(state::TS) where {TS<:AbstractSamplerState}\n\nTrival constraint function\n\n\n\n\n\ntrivial_constraint! - Trivial constraint function \n\nFields\n\nstate - Current state of the sampler\ni     - Index of current iterate\n\n\n\n\n\n","category":"function"},{"location":"sample1/#Recylcing-(EXPERIMENTAL)","page":"Sampling","title":"Recylcing (EXPERIMENTAL)","text":"","category":"section"},{"location":"utils/opts1/#Options-and-Auxiliary-Functions","page":"Options and Auxiliary Functions","title":"Options and Auxiliary Functions","text":"","category":"section"},{"location":"utils/opts1/#Sampler-Options","page":"Options and Auxiliary Functions","title":"Sampler Options","text":"","category":"section"},{"location":"utils/opts1/","page":"Options and Auxiliary Functions","title":"Options and Auxiliary Functions","text":"These options set the number of iterations and the frequency at which data is recorderd.  This is generically used by all samplers.","category":"page"},{"location":"utils/opts1/","page":"Options and Auxiliary Functions","title":"Options and Auxiliary Functions","text":"    MDOptions(; n_iters = 10^4, n_save_iters = 1)","category":"page"},{"location":"utils/opts1/#BasicMD.MDOptions-Tuple{}","page":"Options and Auxiliary Functions","title":"BasicMD.MDOptions","text":"MDOptions(;n_iters = 10^4, n_save_iters=1)\n\nSet options for samplers.\n\nFields\n\nn_iters       - Set the number of iterations of the sampler\nnsaveiters  - Set the frequency at which iterations are saved.  If                 nsaveiters=1, every iteration is saved.  If nsaveiters=n_iters,                 only the final iteration is saved.\n\n\n\n\n\n","category":"method"},{"location":"utils/opts1/#Verlet-Integrator","page":"Options and Auxiliary Functions","title":"Verlet Integrator","text":"","category":"section"},{"location":"utils/opts1/","page":"Options and Auxiliary Functions","title":"Options and Auxiliary Functions","text":"While the main goal of this package is to sample from NVT type ensembles, as it is needed for HMC, the Verlet integrator is included should one wish to sample from the NVE ensemble:","category":"page"},{"location":"utils/opts1/","page":"Options and Auxiliary Functions","title":"Options and Auxiliary Functions","text":"    Verlet","category":"page"},{"location":"utils/opts1/#BasicMD.Verlet","page":"Options and Auxiliary Functions","title":"BasicMD.Verlet","text":"Verlet(∇V!, M, Δt)\n\nSet up the Verlet integrator.\n\nFields\n\n∇V!   - In place gradient of the potential\nM     - Mass (either scalar or vector)\nΔt    - Time step\n\n\n\n\n\n","category":"type"},{"location":"#BasicMD.jl-Documentation","page":"Home","title":"BasicMD.jl Documentation","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a collection of basic routines for Molecular Dynamics simulations written in Julia.  These include","category":"page"},{"location":"","page":"Home","title":"Home","text":"Euler–Maruyama (EM)\nRandom Walk Metropolis (RWM)\nMetropolis Adjusted Langevin (MALA)\nBrünger, Brooks and Karplus (BBK)\nGrønbech-Jensen and Farago (GJF)\nABOBA, BAOAB\nHamiltonian/Hybrid Monte Carlo (HMC)","category":"page"},{"location":"#Caveats","page":"Home","title":"Caveats","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The code assumes that the state space is vector valued.  Thus, even if the problem is one dimensional, you should have initial points and functions formatted appropriately, i.e.","category":"page"},{"location":"","page":"Home","title":"Home","text":"> x0 = [1.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"The mass matrix, M, used in the inertial Langevin integrators and Hamiltonian methods must be diagonal and provided either as a scalar (in the isotropic case) or a vector (in the anisotropic case).  This restriction is in place for performance purposes.\nBBK is currently implemented for a slightly different version of the Langevin SDE than ABOBA/BAOAB.  BBK requires inverting the mass matrix while ABOBA/BAOAB require its square root.\nGJF is implemented in (q,p) coordinates as opposed to (x,v) coordinates. Consequently, the mass term appears slightly differently than in the literature.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Both will sample the associated Boltzmann distribution, but the SDE trajectories will differ when M≂̸I.","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This work was supported in part by the US National Science Foundation Grant DMS-1818716.","category":"page"}]
}