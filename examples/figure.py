"""Module filename: figure.py


"""
import numpy as np
import numpy.fft as nfft


def spectrogramFigure(sigDat, specDat, oName=None):
    r"""Generate spectrogram figure.

    Generate figure showing the intensity normalized spectrogram accompanied by
    subplots that compare the intensities per unit time and frequency, obtained
    for the initial analytic signal, to the normalized marginals, computed for
    the spectrogram.

    Notes
    -----
    Requires the functionality of pythons matplotlib.

    Scales the spectrogram data so that maximum intensity per time and
    frequency is unity.

    In the processing of the input data only the real part of the time domain
    signal Et on input is considered.

    Paramters
    ---------
    sigDat : tuple(array_like, array_like)
        The tuple consists of two arrays in the form (t,Et), where t contains
        the time samples and Et the corrsponding time domain signal.
    res : Results
        The spectrogram ``Results'' data structure with attributes ``tau'' the
        delay time samples for the spectrogram, ``w'' the angular frequency
        samples for the spectrogram, ``P'' the analytic signal spectrogram,
        ``IAE1'' the integrated absolute error of the estimated time marginal,
        and, ``IAE2'' the integrated absolute error of the estimated frequency
        marginal.
    oName : str
        Filename for output figure (default: 'figure0.eps').

    """
    t,Et = sigDat
    tDelay, wOpt, Ptw = specDat.tau, specDat.w, specDat.P

    tMin, tMax = tDelay[0], tDelay[-1]
    wMin, wMax = wOpt[0], wOpt[-1]

    if not oName:
       oName = 'figure0.png'

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.colors as col
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    mm2inch = lambda x: x/10./2.54
    mpl.rcParams['xtick.direction']= 'out'
    mpl.rcParams['ytick.direction']= 'out'
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['axes.labelsize'] = 6
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['figure.figsize'] = mm2inch(90),mm2inch(2./3*90)
    mpl.rcParams["legend.handlelength"] = 1.
    mpl.rcParams["legend.handletextpad"] = 0.15
    mpl.rcParams["legend.borderpad"] = 0.15
    mpl.rcParams["legend.labelspacing"] = 0.15
    cmap=mpl.cm.get_cmap('jet')

    fig = plt.figure()
    plt.subplots_adjust(left=0.12 , bottom=0.17 , right=0.98, top= 0.96, wspace=0.08, hspace=0.06)

    gs = GridSpec(1,1)
    gsSub = GridSpecFromSubplotSpec(4,6,subplot_spec=gs[0,0])
    ax2 = plt.subplot(gsSub[1:,0:5])
    ax1 = plt.subplot(gsSub[0,0:5], sharex=ax2)
    ax3 = plt.subplot(gsSub[1:,5], sharey=ax2)

    def _setColorbar(im, refPos):
        """colorbar helper"""
        x0, y0, w, h = refPos.x0, refPos.y0, refPos.width, refPos.height
        cax = fig.add_axes([x0+0.02*w, y0+0.64*h, 0.02*w, 0.25*h])
        cbar = fig.colorbar(im, cax=cax, ticks=[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],orientation='vertical')
        cbar.outline.set_edgecolor('white')
        cbar.ax.tick_params(color='white', labelcolor='white',right=True, direction='in', labelright=True, labelleft=False, left=False,length=2.)
        cbar.set_ticks((1e-7,1e-5,1e-3,1e-1))

    def _setKey(ax, lns):
        """key helper"""
        labs = [l.get_label() for l in lns]
        lg = ax.legend(lns, labs, title=r"", loc=2, fontsize=6)
        lg = ax.legend_
        lg.get_frame().set_linewidth(0.5)

    def _fft(t, Xt):
        r"""Scaled Fourier transform."""
        return (t[1]-t[0])*nfft.fft(Xt)/np.sqrt(2.*np.pi)

    def _ifft(w, Xw):
        r"""Scaled inverse Fourier transform."""
        return (w[1]-w[0])*nfft.ifft(Xw)/np.sqrt(2.*np.pi)*w.size

    # TOP SUBFIGURE: ##########################################################################
    # COMPARISON OF INTENSITY PER UNIT TIME (ORIGINAL DATA) AND TIME MARGINAL (SPECTROGRAM) ### 
    normEt2 = 1./np.trapz(np.abs(Et)**2,x=t)
    l1 = ax1.semilogy(t,np.abs(Et)**2*normEt2,color='gray',label=r"$|\mathcal{E}(\tau)|^2$")

    P1 = np.trapz(Ptw,x=wOpt,axis=0); normP1 = 1./np.trapz(P1,x=tDelay)
    l2 = ax1.semilogy(tDelay,P1*normP1,color='black',label=r"$\mathsf{P}_1(\tau)$")

    _setKey(ax1,l1+l2)

    ax1.set_ylim(1e-7,2e-2)
    ax1.set_xlim(tMin,tMax)
    ax1.set_yticks((1e-6,1e-4,1e-2))
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.tick_params(axis='y',length=2., direction='in')
    ax1.tick_params(axis='x',length=2., direction='in', labelbottom=False)

    # RIGHT SUBFIGURE: ########################################################################
    # COMPARISON OF INTENSITY PER UNIT FREQUENCY (ORIGINAL DATA) AND FREQUENCY MARGINAL ####### 
    w = nfft.fftfreq(Et.size,d=t[1]-t[0])*2*np.pi
    Ew = _fft(t,np.real(Et))
    Ew[1:Ew.size//2]*=2
    Ew[Ew.size//2+1:]=0
    normEw = 1./np.trapz(np.abs(Ew)**2,x=w)
    l1 = ax3.semilogx(np.abs(Ew)**2*normEw,w,color='gray',label=r"$|\hat{\mathcal{E}}(\omega)|^2$")

    P2 = np.trapz(Ptw,x=tDelay,axis=-1); normP2 = 1./np.trapz(P2,x=wOpt)
    l2 = ax3.semilogx(P2*normP2,wOpt,color='black', label=r"$\mathsf{P}_2(\omega)$")

    _setKey(ax3,l1+l2)

    ax3.set_ylim(wMin,wMax)
    ax3.set_xlim(5e-5,20)
    ax3.set_xticks((1e-4,1e-2,1e0))
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.tick_params(axis='y',length=2., direction='in', labelleft=False)
    ax3.tick_params(axis='x',length=2., direction='in')

    # CENTER FIGURE: SPECTROGRAM DATA ##########################################################
    I = Ptw[:-1,:-1]/Ptw.max()
    imA = ax2.pcolorfast(tDelay,wOpt,I,norm=col.LogNorm(vmin=1e-5*I.max(),vmax=I.max()),cmap=cmap)

    _setColorbar(imA,ax2.get_position())

    ax2.set_xlim(tMin,tMax)
    ax2.set_ylim(wMin,wMax)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.tick_params(axis='y',length=2., direction='in', color='white')
    ax2.tick_params(axis='x',length=2., direction='in', color='white')
    ax2.set_xlabel(r"Delay $\tau\,\mathrm{(fs)}$")
    ax2.set_ylabel(r"Angular frequency $\omega\,\mathrm{(rad/fs)}$")

    ax2.text(0.02,0.975, r"$P_S(\tau,\omega)$",color='white',
      horizontalalignment='left',
      verticalalignment='top',
      transform=ax2.transAxes,fontweight='heavy')

    plt.savefig(oName,dpi=800)
    plt.close()


# EOF: figure.py
