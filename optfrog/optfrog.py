"""Module filename: optfrog.py

Implements algorithms and data structures for the computation of analytic
signal spectrograms with opimized time and frequency resolution.

"""
import sys
import numpy as np
import numpy.fft as nfft
import scipy.optimize as so


class Results(object):
    r"""Represents the spectrogram data."""
    pass

def _fft(t, Xt):
    r"""Scaled Fourier transform."""
    return (t[1]-t[0])*nfft.fft(Xt)/np.sqrt(2.*np.pi)

def _ifft(w, Xw):
    r"""Scaled inverse Fourier transform."""
    return (w[1]-w[0])*nfft.ifft(Xw)/np.sqrt(2.*np.pi)*w.size

def _analyticSignal(t, Et):
    r"""Construction of analytic signal for input field."""
    w = nfft.fftfreq(t.size,d=t[1]-t[0])*2*np.pi
    Zw = _fft(t,np.real(Et))
    Zw[1:Zw.size//2]*=2
    Zw[Zw.size//2+1:]=0
    Zt = _ifft(w,Zw)
    return (t,Zt),(w,Zw)

def _L2Norm(t, Et):
    r"""L2 norm for input field."""
    return np.trapz(np.abs(Et)**2,x=t)


def _zeroPad(q):
    r"""Zero padding for input array."""
    return np.pad(q,(0,2**(q.size-1).bit_length()-q.size),mode='constant')


def vanillaFrog(t, Et, hFunc, tLim = None, wLim = None ):
    r"""Spectrogram for time-domain analytic signal.

    Compute an analytic signal spectrogram for time-domain analytic signal. The
    resulting time-frequency representation allows to interpret the
    time-frequency characteristics of the underlying signal [1]_.

    Parameters
    ----------
    t : array_like
        Array containing time samples.
    Et : array_like
        Array containing analytic signal in time-domain.
    hFunc : object
        Custom window function for short-time fourier transform.
    tLim : tuple(float, float, int)
        Boundary conditions for filtering and slicing of time samples in the
        form (tMin, tMax, mt), specifying lower bound for time samples, upper
        bound for time samples and an integer for filtering every mt-th time
        sample.
    wLim : tuple(float, float, int)
        Boundary conditions for filtering and slicing of time samples in the
        form (wMin, wMax, mw), specifying lower bound for frequency samples,
        upper bound for frequency samples and an integer for filtering every
        mw-th frequency sample.

    Returns
    -------
    res : Results
        The spectrogram ``Results'' data structure with attributes ``tau'' the
        delay time samples for the spectrogram, ``w'' the angular frequency
        samples for the spectrogram, ``P'' the analytic signal spectrogram,
        ``IAE1'' the integrated absolute error of the estimated time marginal,
        and, ``IAE2'' the integrated absolute error of the estimated frequency
        marginal.

    Notes
    -----
    Analytic signal for the real part of the input field is obtained using
    their respective frequency domain definitions [2]_.

    The spectrogram is obtained as the squared magnitude of a short-time
    fourier-transform (STFT) using a custom window function [3]_.

    STFT resolution issues: narrow window-function yields good time resolution
    but poor frequency resolution; wide window function yields good frequency
    resolution but bad time resolution.

    Time-domain signal is zero padded so as to allow for efficient processing
    via FFT.

    Synonymous with frequency resolved electrical gating (FREG) analysis [4]_
    and cross-correlation frequency resolved optical gating (XFROG) analysis
    [5]_.

    In the field of nonlinear optics the respective analysis is referred to as
    frequency resolved optical gating (i.e. frog) analysis. The routine
    vanillaFrog computes a straight forward spectrogram for a particular window
    function, hence the preceding term vanilla.

    References
    ----------
    .. [1] Cohen, L. Time-Frequency distributions - A Review. Proceedings of
       the IEEE, 77 (1989) 941.

    .. [2] Marple, S L. Computing the discrete-time 'analytic' signal via FFT.
       IEEE Trans. Signal Processing, 47 (1999) 2600.

    .. [3] https://en.wikipedia.org/wiki/Short-time_Fourier_transform

    .. [4] Dorrer, C and Kang, I. Simultaneous temporal characterization
       of telecommunication optical pulses and modulators by use of spectrograms
       Opt. Lett., 27 (2002) 1315.

    .. [5] Linden, S and Giessen, H and Kuhl, J. XFROG - A New Method for
       Amplitude and Phase Characterization of Weak Ultrashort Pulses. Phys. Stat.
       Sol. 206 (1998) 119.

    """
    if tLim is None:
        tLim = (-np.inf, np.inf, 1)
    if wLim is None:
        wLim = (0., np.inf, 1)
    # INITIALIZE RESULTS ################################################################
    res = Results()
    # ZERO PAD TO FACILITATE EFFICIENT USE OF FFT #######################################
    t = _zeroPad(t)
    Et = _zeroPad(Et)
    # NORMALIZE INPUT FIELD SO THAT L2-NORM OF AS IS UNITY ##############################
    Et = np.real(Et)/np.sqrt(2.*_L2Norm(t,np.real(Et)))
    # OBTAIN ANALYTIC SIGNAL USING ITS FREQUENCY DOMAIN FORMULATION ##################### 
    (t,Zt),(w,Zw) = _analyticSignal(t,Et)
    w = nfft.fftshift(w)
    Zw = nfft.fftshift(Zw)

    # OBTAIN DOWNSAMPLED VARIABLE RANGES ################################################
    tMin,tMax,mt = tLim
    tMask = np.logical_and(t>tMin,t<tMax)
    tau = t[tMask]; tau=tau[::mt]
    Zts = Zt[tMask]; Zts=Zts[::mt]

    wMin,wMax,mw = wLim
    wMask = np.logical_and(w>wMin,w<wMax)
    wOpt = w[wMask]; wOpt=wOpt[::mw]
    Zws = Zw[wMask]; Zws=Zws[::mw]

    # INITIALIZE SPECTROGRAM ############################################################
    P_STFT = np.zeros((wOpt.size, tau.size))
    # INTERATE OVER DOWNSAMPLED DELAY TIMES ############################################# 
    for idx,tauCurr in enumerate(tau):
        # COMPUTE SHORT-TIME FOURIER TRANSFORM OF WINDOWED SIGNAL #######################
        STFT = nfft.fftshift(_fft(t,Zt*hFunc(t-tauCurr)))
        # OBTAIN BOUNDED AND DOWNSAMPLED SPECTROGRAM FROM STFT ##########################  
        P_STFT[:,idx] = np.abs(STFT[wMask][::mw])**2

    # UPDATE OUTPUT DATA STRUCTURE WITH PRELIMINARY RESULTS #############################
    res.tau = tau
    res.w = wOpt
    res.P = P_STFT

    # COMPUTE SUMMARY MEASURES TO ASSESS SPECTROGRAM FIT ################################
    P1 = timeMarginal(res)
    P2 = frequencyMarginal(res)
    Es = totalEnergy(res)

    # QUANTIFY FIT OF SPECTROGRAM TO ORIGINAL DATA ######################################
    res.IAE1 = np.trapz(np.abs(np.abs(Zts)**2-P1/Es),x=tau)
    res.IAE2 = np.trapz(np.abs(np.abs(Zws)**2-P2/Es),x=wOpt)
    return res


def timeMarginal(res):
    r"""Spectrogram based time marginal.

    Compute marginal distribution in time based on the spectrogram.

    Parameters
    ----------
    res : Results
        The spectrogram ``Results'' data structure with attributes ``tau'' the
        delay time samples for the spectrogram, ``w'' the angular frequency
        samples for the spectrogram, ``P'' the analytic signal spectrogram,
        ``IAE1'' the integrated absolute error of the estimated time marginal,
        and, ``IAE2'' the integrated absolute error of the estimated frequency
        marginal.

    Returns
    -------
    P1 : array_like
        Time marginal estimated from the spectrogram.

    References
    ----------
    .. [1] Cohen, L. Time-Frequency distributions - A Review. Proceedings of
       the IEEE, 77 (1989) 941.

    """
    return np.trapz(res.P,x=res.w,axis=0)


def frequencyMarginal(res):
    r"""Spectrogram based frequency marginal.

    Compute marginal distribution in frequency based on the spectrogram.

    Parameters
    ----------
    res : Results
        The spectrogram ``Results'' data structure with attributes ``tau'' the
        delay time samples for the spectrogram, ``w'' the angular frequency
        samples for the spectrogram, ``P'' the analytic signal spectrogram,
        ``IAE1'' the integrated absolute error of the estimated time marginal,
        and, ``IAE2'' the integrated absolute error of the estimated frequency
        marginal.

    Returns
    -------
    P2 : array_like
        Frequency marginal estimated from the spectrogram.

    References
    ----------
    .. [1] Cohen, L. Time-Frequency distributions - A Review. Proceedings of
       the IEEE, 77 (1989) 941.

    """
    return np.trapz(res.P,x=res.tau,axis=-1)


def totalEnergy(res):
    r"""Total energy.

    Compute total energy provided by the spectrogram approximation of the time-frequency
    characteristics of the signal.

    Parameters
    ----------
    res : Results
        The spectrogram ``Results'' data structure with attributes ``tau'' the
        delay time samples for the spectrogram, ``w'' the angular frequency
        samples for the spectrogram, ``P'' the analytic signal spectrogram,
        ``IAE1'' the integrated absolute error of the estimated time marginal,
        and, ``IAE2'' the integrated absolute error of the estimated frequency
        marginal.

    Returns
    -------
    Es : float
        Total energy estimated from the spectrogram.

    Notes
    -----
    Since a spectrogram does not formally satisfy the marginals, the total
    energy provided by the spectrogram generally differs from the total energy
    of the original signal.

    References
    ----------
    .. [1] Cohen, L. Time-Frequency distributions - A Review. Proceedings of
       the IEEE, 77 (1989) 941.

    """
    return np.trapz(timeMarginal(res),x=res.tau)


def optFrog(t,Et, hFunc, alpha=0., tLim=None, wLim=None, disp=True):
    r"""Time-frequency relolution optimized spectrogram for time-domain
    analytic signal.

    Compute an optimized analytic signal spectrogram for time-domain analytic
    signal. The resulting time-frequency representation allows to interpret the
    time-frequency characteristics of the underlying signal [1]_.  The obtained
    spectrograms are optimal in the sense that they use an optimal parameter
    setting for the user supplied window function that minimizes the integrated
    absolute error between the intensities per unit time and frequancy,
    obtained from the initial analytic signal, and the normalized marginals,
    computed from the spectrogram.

    Parameters
    ----------
    t : array_like
        Array containing time samples.
    Et : array_like
        Array containing analytic signal in time-domain.
    hFunc : object
        Custom window function for short-time fourier transform.
    alpha : float
        Parameter allowing to adjust relative importance of time or frequency
        resolution (default alpha=0.0).
    tLim : tuple(float, float, int)
        Boundary conditions for filtering and slicing of time samples in the
        form (tMin, tMax, mt), specifying lower bound for time samples, upper
        bound for time samples and an integer for filtering every mt-th time
        sample.
    wLim : tuple(float, float, int)
        Boundary conditions for filtering and slicing of time samples in the
        form (wMin, wMax, mw), specifying lower bound for frequency samples,
        upper bound for frequency samples and an integer for filtering every
        mw-th frequency sample.
    disp : bool
        Set to True to display convergece info (default: True).

    Returns
    -------
    res : Results
        The spectrogram ``Results'' data structure with attributes ``tau'' the
        delay time samples for the spectrogram, ``w'' the angular frequency
        samples for the spectrogram, ``P'' the analytic signal spectrogram,
        ``IAE1'' the integrated absolute error of the estimated time marginal,
        and, ``IAE2'' the integrated absolute error of the estimated frequency
        marginal.

    Notes
    -----
    The optimization procedure uses the function vanillaFrog() for the
    calculation of the spectrograms.

    As discussed in Ref. [6]_, if a window function approximates one marginal
    closely but the other one poorly, the spectrogram might appear distorted
    and unreasonable. In contranst, visual inspection confirms that
    spectrograms with the above property yield a reasonable time-frequency
    characterization of the signal under scrutiny.

    References
    ----------
    .. [1] Cohen, L. Time-Frequency distributions - A Review. Proceedings of
       the IEEE, 77 (1989) 941.

    .. [2] Marple, S L. Computing the discrete-time 'analytic' signal via FFT.
       IEEE Trans. Signal Processing, 47 (1999) 2600.

    .. [3] https://en.wikipedia.org/wiki/Short-time_Fourier_transform

    .. [4] Dorrer, C and Kang, I. Simultaneous temporal characterization
       of telecommunication optical pulses and modulators by use of spectrograms
       Opt. Lett., 27 (2002) 1315.

    .. [5] Linden, S and Giessen, H and Kuhl, J. XFROG - A New Method for
       Amplitude and Phase Characterization of Weak Ultrashort Pulses. Phys. Stat.
       Sol. 206 (1998) 119.

    .. [6] Cohen, L and Loughlin, P J. The marginals and time-frequency
       distributions. Proc. SPIE 5205, Advanced Signal Processing Algorithms,
       Architechtures, and Implementations XIII, doi: 10.1117/12.51389.

    """
    if tLim is None:
        tLim = (-np.inf, np.inf, 1)
    if wLim is None:
        wLim = (0., np.inf, 1)
    # SET REASONABLE BOUNDS FOR WIDTH PARAMETER OF WINDOW FUNCION #######################
    s0Min = 10.*(t[1]-t[0])
    s0Max = 0.25*(min([t.max(),tLim[1]])-max([t.min(),tLim[0]]))

    # NORMALIZE INPUT FIELD SO THAT L2-NORM OF AS IS UNITY ##############################
    Et = np.real(Et)/np.sqrt(2.*_L2Norm(t,np.real(Et)))

    if disp:
      sys.stderr.write("# optFrog convergence info:\n")
      sys.stderr.write("# (sigma) (Q(sigma, alpha=%lf)) (IAE1) (IAE2)\n"%(alpha))

    def _objFunc(s0,alpha, t, Et, tLim, wLim):
        res = vanillaFrog(t, Et, hFunc(s0), tLim, wLim)
        QVal = 0.5*(1.+alpha)*res.IAE1 + 0.5*(1.-alpha)*res.IAE2
        if disp:
          sys.stderr.write("%lf %lf %lf %lf\n"%(s0, QVal, res.IAE1, res.IAE2))
        return QVal

    optRes = so.minimize_scalar( lambda x: _objFunc(x, alpha, t, Et, tLim, wLim), bounds=(s0Min,s0Max), method='bounded')

    # YIELD SPECTROGRAM WITH OPTIMAL WIDTH PARAMETER ####################################
    return vanillaFrog(t, Et, hFunc(optRes.x), tLim, wLim)


# EOF: optfrog.py
