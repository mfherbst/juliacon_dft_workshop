# Solution to the exercise "Computing stresses"
# in the notebook 6_Floating_point_error.ipynb

using DFTK
using GenericLinearAlgebra  # Needed to enable generic fallbacks in DFTK
using DoubleFloats          # Defines Double64

# Silicon lattice
a = 10.26
lattice = a / 2 .* [[0 1 1.]; [1 0 1.]; [1 1 0.]]
Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
atoms = [Si => [ones(3)/8, -ones(3)/8]]

# Repeat the Float32 calculation here for the sake of having everything
# in one place.
println("Running Float32 ...")
model = model_DFT(Array{Float32}(lattice), atoms, [:lda_x, :lda_c_vwn])
basis = PlaneWaveBasis(model; Ecut=7, kgrid=[1, 1, 1], fft_size=(16, 16, 16))
scfres = self_consistent_field(basis, tol=1e-16, maxiter=40)
results = Dict{DataType, Any}(Float32 => scfres.energies.total)

println("\nRunning Float64 ...")
model = model_DFT(Array{Float64}(lattice), atoms, [:lda_x, :lda_c_vwn])
basis = PlaneWaveBasis(model; Ecut=7, kgrid=[1, 1, 1], fft_size=(16, 16, 16))
scfres = self_consistent_field(basis, tol=1e-16, maxiter=40)
results[Float64] = scfres.energies.total


println("\nRunning Double64 ...")
model = model_DFT(Array{Double64}(lattice), atoms, [:lda_x, :lda_c_vwn])
basis = PlaneWaveBasis(model; Ecut=7, kgrid=[1, 1, 1], fft_size=(16, 16, 16))
scfres = self_consistent_field(basis, tol=1e-16, maxiter=40)
results[Double64] = scfres.energies.total

println()
println("Float32:   $(results[Float32])")
println("Float64:   $(results[Float64])")
println("Double64:  $(results[Double64])")

println()
println("Errors versus Double64:")
println("Float32:   $(results[Float32] - results[Double64])")   # Order of 1e-7
println("Float64:   $(results[Float64] - results[Double64])")   # Order of 1e-14