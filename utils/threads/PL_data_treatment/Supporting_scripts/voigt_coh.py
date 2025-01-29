import numpy as np

def Voigt_coh(x,A0,A1,A2,A3,A4,A5,A6):

    # Parameters:
    # X: delay
    # A[0]: Amplitude
    # A[1]: Collision time T2 related to Gaussian broadening
    # A[2]: Position corresponding to maximum visibility(0 delay)
    # A[3]: Coherence time T2
    #
    # Optional
    # parameters:
    # A[4]: Offset – should be set to 0
    # A[5]: Linear term
    # A[6]: Quadratic term
    #
    # Derived
    # parameters(optional):
    # Linewidht_Gauss: 2 * hbar / (A[1]) * sqrt(2 * alog(2.0));
    # FWHM in ueV(Gauss):
    # Linewidht_Lorentz: 2 * hbar / A[3];
    # FWHM in ueV(Lorentz)
    # Linewidth(ueV) = dEl * c1 + sqrt(dEl ^ 2 * c2 + dEg ^ 2);
    # FWHM in ueV(Voigt)

    return A0 * np.exp(-1 / 2 * ((x - A2) / A1) ** 2) * np.exp(-1 * np.abs(x - A2) / A3) + A4 + A5*x + A6*x**2