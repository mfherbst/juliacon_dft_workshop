# Solution to the exercise computing the Kerker-preconditioned eigenvalues
# of aluminium in 5_Analysing_SCF_convergence.jl

# Copy to make script self-contained
using DFTK
using LinearAlgebra
using Plots
using Statistics
setup_threading()

function aluminium_setup(repeat=1; Ecut=7.0, kgrid=[1, 1, 1])
    a = 7.65339
    lattice = diagm(fill(a, 3))
    Al = ElementPsp(:Al, psp=load_psp("hgh/lda/al-q3"))
    atoms = [Al => [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]]

    # Make supercell in pymatgen
    mg_struct = pymatgen_structure(lattice, atoms)
    mg_struct.make_supercell([1, 1, repeat])
    lattice = load_lattice(mg_struct)
    atoms = [Al => [s.frac_coords for s in mg_struct.sites]];

    # Construct the model
    model = model_LDA(lattice, atoms, temperature=1e-3, symmetries=false)
    PlaneWaveBasis(model; Ecut, kgrid)
end

function construct_Pinv_epsilon(scfres)
    basis = scfres.basis
    
    Pinv_Kerker(δρ) = DFTK.mix_density(KerkerMixing(), basis, δρ)

    function epsilon(δρ)  # Apply ε† = 1 - χ0 Khxc
        δV   = apply_kernel(basis, δρ; ρ=scfres.ρ)
        χ0δV = apply_χ0(scfres.ham, scfres.ψ, scfres.εF, scfres.eigenvalues, δV)
        δρ - χ0δV   
    end    
    
    epsilon, Pinv_Kerker
end

function plot_mode(mode)
    # Average along x axis
    mode_yz = mean(real.(mode), dims=1)[1, :, :, 1]
    heatmap(mode_yz, c=:RdBu_11, aspect_ratio=1, grid=false,
            legend=false, clim=(-0.006, 0.006))
end

# Start of new stuff ...
using KrylovKit

repeat = 3
scfres = self_consistent_field(aluminium_setup(repeat); tol=1e-12, mixing=KerkerMixing());
epsilon, Pinv_Kerker = construct_Pinv_epsilon(scfres)
operator = Pinv_Kerker ∘ epsilon   # ∘ = \circ<Tab>

λ_large, X_large, info = eigsolve(operator, randn(size(scfres.ρ)), 4, :LM;
                                  tol=1e-4, eager=true, verbosity=2)
@assert info.converged ≥ 4
λ_max = maximum(real.(λ_large))
println("Largest eigenvalue: $(λ_max)")

# Getting the smallest eigenvalue here is quite a pain ... we will just assume it is 0.8
# λ_small, X_small, info = eigsolve(epsilon, randn(size(scfres.ρ)), 2, EigSorter(abs, rev=false);
#                                   tol=1e-3, eager=true, verbosity=2)
# @assert info.converged ≥ 2
# λ_min = minimum(real.(λ_small))
λ_min = 0.8
cond = λ_max / λ_min

println("Smallest eigenvalue: $(λ_min)")  # about 0.8
println("Largest eigenvalue:  $(λ_max)")  # about 8.6
println("Condition number:    $cond")     # about 10

# Plot largest mode to pdf
savefig(plot_mode(X_large[1]), "5_ex1_aluminium_Kerker_repeat_$repeat.pdf")

# So by using the Kerker preconditioner for the repeat 3 the condition number reduced from 30 to 10.
#
# For repeat 6:
#     largest eigenvalue without Kerker:  100
#     largest eigenvalue with    Kerker:  5.5
# So by using Kerker the condition number is reduced from 105 to 6.8.