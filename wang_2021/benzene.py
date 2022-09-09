INITIAL_GEOMETRY = """
C   1.209     0.698   0.000
C   0.000     1.396   0.000
C  -1.209     0.698   0.000
C  -1.209    -0.698   0.000
C  -0.000    -1.396   0.000
C   1.209    -0.698   0.000
H   2.159     1.247   0.000
H   0.000     2.493   0.000
H  -2.159     1.247   0.000
H  -2.159    -1.247   0.000
H  -0.000    -2.493   0.000
H   2.159    -1.247   0.000
""" #in Angstrom
WANG_XC = "LDA"


def dipole_moments(initial_geometry:str,xc_functional:str,number_of_states:int):
    """returns list of dipole moments"""
    from pyscf import gto, dft

    mol = gto.M(atom = initial_geometry,
                basis='6-31g',
                symmetry=True)
    mf = dft.RKS(mol)
    mf.xc = xc_functional

    from pyscf.geomopt.geometric_solver import optimize
    conv_params = {
        'convergence_energy': 1e-6,  # Eh
        'convergence_grms': 3e-4,    # Eh/Bohr
        'convergence_gmax': 4.5e-4,  # Eh/Bohr
        'convergence_drms': 1.2e-3,  # Angstrom
        'convergence_dmax': 1.8e-3,  # Angstrom
    }
    mol_eq = optimize(mf, **conv_params)

    mf = dft.RKS(mol_eq)
    mf.xc = xc_functional
    mf.kernel()

    from pyscf import tddft #scf
    mytd = tddft.TDDFT(mf)
    mytd.nstates = number_of_states
    energies = mytd.kernel()
    mytd.analyze(verbose=4)
    return energies,mytd.transition_dipole()