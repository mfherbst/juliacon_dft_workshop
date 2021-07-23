# Solution to the exercise computing the Helium dielectric eigenvalues
# in 5_Analysing_SCF_convergence.jl

# Copy to make script self-contained
using DFTK
using LinearAlgebra
using KrylovKit

function helium_setup(repeat=40; Ecut=7.0, kgrid=[1, 1, 1])
    a = 5
    lattice = diagm(fill(a, 3))
    He = ElementPsp(:He, psp=load_psp("hgh/lda/he-q2"))
    atoms = [He => [[0.0, 0.0, 0.0]]]

    # Make supercell in pymatgen
    mg_struct = pymatgen_structure(lattice, atoms)
    mg_struct.make_supercell([1, 1, repeat])
    lattice = load_lattice(mg_struct)
    atoms = [He => [s.frac_coords for s in mg_struct.sites]];

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

scfres = self_consistent_field(helium_setup(30); tol=1e-12);
epsilon, Pinv_Kerker = construct_Pinv_epsilon(scfres)

for (label, operator) in [("ε", epsilon), ("P⁻¹ε", Pinv_Kerker ∘ epsilon)]
    println()
    println("--------")
    println()
    println("Compute λ_large for $label")
    λ_large, X_large, info = eigsolve(operator, randn(size(scfres.ρ)), 2, :LM;
                                      tol=1e-3, eager=true, verbosity=2)
    @assert info.converged ≥ 2
    λ_max = maximum(real.(λ_large))
        
    println("Compute λ_small for $label")
    λ_small, X_small, info = eigsolve(operator, randn(size(scfres.ρ)), 2, EigSorter(abs, rev=false);
                                      tol=1e-3, eager=true, verbosity=2)
    @assert info.converged ≥ 2
    λ_min = minimum(real.(λ_small))
        
    cond = λ_max / λ_min
    println()
    println("Smallest eigenvalue $label: $(λ_min)")
    println("Largest eigenvalue $label:  $(λ_max)")
    println("Condition number $label:    $cond")
end


# Example output:
#
# Smallest eigenvalue ε: 0.9517149905415806
# Largest eigenvalue ε:  1.4438934058378814
# Condition number ε:    1.517148957605704
#
#
# Smallest eigenvalue P⁻¹ε: 0.0036831161341721696
# Largest eigenvalue P⁻¹ε:  1.1491558189981426
# Condition number P⁻¹ε:    312.0064035820665
#
# Using the Kerker mixing on the Helium system has little effect on the largest eigenvalue, but it
# dramatically lowers the smallest eigenvalue (from 0.95 to about 3e-2) and as a result the condition
# number *increases* from a dream value of 1.5 to about 300. As a result SCF convergence is much worse
# when using the Kerker mixing on this insulating system.