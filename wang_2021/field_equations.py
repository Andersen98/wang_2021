
import numpy as np
def field_strength_of_mode(k, lambda_c, number:int):
    """
    Inputs:
        k (3d vector): desired mode
        lambda_c (3d vector): center mode strength 
        wc (float): center freqency of mode 
        n(int): number of modes 
    Returns:
        field_strength (float): Field strength at a mode k given Lorenzian 
                        distribution about wc with magnitude lambda_c.
                        The sum over the n equally spaced modes adds 
                        up to lambda_c.
    Sources: 
        Equations  3 and 4 of
        `Wang, Derek S., Tomáš Neuman, Johannes Flick, and Prineha Narang.
         “Light-Matter Interaction of a Molecule in a Dissipative Cavity 
         from First Principles.” Journal of Chemical Physics 154, no. 10 (2021). 
         https://doi.org/10.1063/5.0036283.`
    """
    pass

def field_strength_of_mode(energy,volume,epsilon):
    """
    Gives the field strength magnitude of a field mode 
    given energy and volume.

    **NOTE THAT THIS IS 3D CASE**

    Inputs:
        energy (float): positive real, can be hbar omega, or h*f
        volume (float): positive real, 3d volume of mode
        epsilon (float): whatever epsilon is in your units/system

    Outputs:
        field magnitude (float): Equation 1.1.20 of Scully and Zubairy, Quantum Optics.
    """
    return np.sqrt(energy/(2*epsilon*volume))