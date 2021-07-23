# Solution to the exercise "Computing stresses"
# in the notebook 7_Properties_automatic_differentation.ipynb

using ForwardDiff
using DFTK
scfres = compute_silicon(zeros(3))

function recompute_silicon_energy_stresses(lattice)
    atoms = scfres.basis.model.atoms
    new_model  = model_DFT(lattice, atoms, [:lda_x, :lda_c_vwn]; symmetries=false)
    new_basis  = PlaneWaveBasis(new_model; Ecut=13, kgrid=[2, 2, 2]);
    ρ = DFTK.compute_density(new_basis, scfres.ψ, scfres.occupation)
    energy_hamiltonian(new_basis, scfres.ψ, scfres.occupation; ρ=scfres.ρ).E.total
end

L = scfres.basis.model.lattice

# Just tranform the stress definition into Julia code:
ForwardDiff.gradient(M -> recompute_silicon_energy_stresses((I + M) * L), zero(L)) / det(L)