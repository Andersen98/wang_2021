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

def cavity_strength(volume):
    """
    Cavity strength as a function of volume.
    Solves eqns 3,4,
    """
