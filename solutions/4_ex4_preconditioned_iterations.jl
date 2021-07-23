# Solution for the exercise "Preconditioned iterations"
# in 4_Solving_the_SCF_problem.ipynb

# Setup copied to have a stand-alone script
using DFTK
using LinearAlgebra
setup_threading()

function aluminium_setup(repeat=1; Ecut=13.0, kgrid=[2, 2, 2])
    a = 7.65339
    lattice = diagm(fill(a, 3))
    Al = ElementPsp(:Al, psp=load_psp("hgh/lda/al-q3"))
    atoms = [Al => [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]]

    # Make supercell in pymatgen:
    # We convert our lattice to the conventions used in pymatgen
    # and then back ...
    mg_struct = pymatgen_structure(lattice, atoms)
    mg_struct.make_supercell([1, 1, repeat])
    lattice = load_lattice(mg_struct)
    atoms = [Al => [s.frac_coords for s in mg_struct.sites]];

    # Construct the LDA model and discretise
    model = model_LDA(lattice, atoms, temperature=1e-3)
    PlaneWaveBasis(model; Ecut, kgrid)
end

function run_scf(repeat, mixing)
    scfres = self_consistent_field(aluminium_setup(repeat);
                                   damping=0.8, mixing=mixing(), tol=1e-6)
    scfres.n_iter
end

using DataFrames
df = DataFrame(mixings=Symbol[], repeats=Int[], n_iter=Int[])
for repeat in [1, 2, 4, 6, 8]
    for mixing in [KerkerMixing, DielectricMixing, LdosMixing]
        println("$mixing  --- $repeat")
        push!(df, (Symbol(mixing), repeat, run_scf(repeat, mixing)))
        println()
        println()
    end
end

println(df)
#15×3 DataFrame
# Row │ mixings           repeats  n_iter 
#     │ Symbol            Int64    Int64  
#─────┼───────────────────────────────────
#   1 │ KerkerMixing            1       5
#   2 │ DielectricMixing        1       5
#   3 │ LdosMixing              1      11
#
#   4 │ KerkerMixing            2       5
#   5 │ DielectricMixing        2       5
#   6 │ LdosMixing              2       7
#
#   7 │ KerkerMixing            4       6
#   8 │ DielectricMixing        4       7
#   9 │ LdosMixing              4       9
#
#  10 │ KerkerMixing            6       8
#  11 │ DielectricMixing        6      12
#  12 │ LdosMixing              6       9
#
#  13 │ KerkerMixing            8       8
#  14 │ DielectricMixing        8      12
#  15 │ LdosMixing              8       9

# For the small system sizes LdosMixing does a little worse, but between Dielectric
# and Kerker there is basically no difference, since for the small sizes the system-specific
# details of aluminium do not yet become visible.

# As the system size increases, LdosMixing catches up and becomes as good as Kerker.
# Dielectric (which is best suited for semiconductors) becomes less and less able to help
# the convergence and the number of required iterations increases slowly with system size.
