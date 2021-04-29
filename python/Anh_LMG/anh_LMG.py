def Expected_n(Sval, avecs_U1):
    '''
    
    Expected value of n for eigenvectors in the U(1) dynamical symmetry basis
    
    '''
    import numpy as np
    n_vector = Sval + np.arange(-Sval,Sval+1)
    return np.sum(((avecs_U1**2).T*n_vector).T, axis=0)
################################################################################
def Expected_n_parity(Sval, avecs_U1_parity, parity = "Plus"):
    '''
    
    Expected value of n for well defined parity eigenvectors in the U(1) dynamical symmetry basis
    
    '''
    import numpy as np
    # Only even Sval values
    if parity == "Plus":
        n_vector = Sval + np.arange(-Sval,Sval+1,2)
    elif parity == "Minus":
        n_vector = Sval + np.arange(-Sval+1,Sval,2)

    return np.sum(((avecs_U1_parity**2).T*n_vector).T, axis=0)
################################################################################
def IPR(avecs_U1_parity):
    '''
    
    Inverse participation ration for eigenvectors in the U(1) dynamical symmetry basis
    
    '''
    import numpy as np
    #
    return 1/np.sum(avecs_U1_parity**4, axis=0)
################################################################################
