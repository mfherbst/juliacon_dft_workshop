# Solution for the exercise "Damped iterations"
# in 4_Solving_the_SCF_problem.ipynb

# Setup copied to have a stand-alone script
using DFTK
using LinearAlgebra

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


function damped_fixed_point_iteration(F, ρ₀, maxiter; tol)
    # F:        The SCF step function
    # ρ₀:       The initial guess density
    # maxiter:  The maximal number of iterations to be performed
    # tol:      The selected convergence tolerance
    
    α = 0.77  # <-- Set damping value here
    
    ρ  = ρ₀
    Fρ = F(ρ)
    for n = 1:maxiter 
        # If change less than tolerance, break iterations:
        if norm(Fρ - ρ) < tol
            break
        end

        Rρ    = Fρ - ρ      # Compute residual
        ρnew  = ρ + α * Rρ  # Damped update
        
        # Compute next step
        ρ  = ρnew
        Fρ = F(ρ)
    end
    
    # Return some stuff DFTK needs ...
    (fixpoint=ρ, converged=norm(Fρ - ρ) < tol)
end;

# use this algorithm inside DFTK's SCF
self_consistent_field(aluminium_setup();
                      solver=damped_fixed_point_iteration,
                      damping=1.0,
                      maxiter=40)
nothing


# Very fast convergence is obtained for about α = 0.75 in only 6 SCF steps
# Already slightly larger values break convergence (e.g. α = 0.78 fails to converge)