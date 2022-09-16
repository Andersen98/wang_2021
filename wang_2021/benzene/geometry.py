import numpy as np
from scipy.spatial.transform import Rotation as R
from functools import reduce
BENZENE_CC_ANGSTROM = 1.39
BENZENE_CH_ANGSTROM = 1.09

def initial_structure(cc_angstrom:float,ch_angstrom:float) -> str:
    """Gives bezene structure in xyz format.

    Args:
        cc_angstrom: Carbon-Carbon bond length in angstroms.
        ch_angstrom: Carbon-Hydrogen bond length in angstroms.

    Returns:
        Benzene xyz string.
    """
    r_carbon = cc_angstrom
    r_hydrogen = r_carbon + ch_angstrom

    num_rots=6 #6 fold rotational symmetry
    delta_theta = 2*np.pi/6
    rotations = [R.from_euler('z',n*delta_theta) for n in range(num_rots)]
    carbons = [r.apply(r_carbon*np.array([1.0,0.0,0.0])) for r in rotations]
    hydrogrens= [r.apply(r_hydrogen*np.array([1.0,0.0,0.0])) for r in rotations]

    carbon_strs=['C {x:.2f} {y:.2f} {z:.2f}'.format(x=vec[0],y=vec[1],z=vec[2]) for vec in carbons]
    hydrogen_strs=['H {x:.2f} {y:.2f} {z:.2f}'.format(x=vec[0],y=vec[1],z=vec[2]) for vec in hydrogrens]
    benzene_geometry = carbon_strs+ hydrogen_strs
    benzene_geometry_str = reduce(lambda x,y:x+'\n'+y, benzene_geometry)
    return benzene_geometry_str
