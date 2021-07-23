# Solution to the exercise "Modelling aluminium"
# in the notebook 3_Density_functional_theory.ipynb

using DFTK
using LinearAlgebra

a = 7.65339
lattice = diagm(fill(a, 3))
Al = ElementPsp(:Al, psp=load_psp("hgh/lda/al-q3"))
atoms = [Al => [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]]

# Without smearing (won't work!)
# model = model_LDA(lattice, atoms)

# With smearing: Setup model and discretise
model = model_LDA(lattice, atoms, temperature=1e-3)
basis = PlaneWaveBasis(model; Ecut=15, kgrid=[2, 2, 2])

# Run the SCF
scfres = self_consistent_field(basis)

# Plot the bands
plot_bandstructure(scfres)