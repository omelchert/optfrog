"""Script filename: example_optFrog.py

Exemplary calculation of an optFrog trace for data obtained from
the numerical propagation of a short and intense few-cycle optical
pulse in presence of the refractive index profile of an endlessly single
mode photonic crystal fiber.

"""
import sys
import numpy as np
import numpy.fft as nfft
from optfrog import optFrog
from figure import spectrogramFigure


def main():

    tMin= -500.0
    tMax= 5800.0
    wMin= 0.75
    wMax= 3.25
    s0 = 36.0
    a0 = 0.0
    fName = './data/exampleData_pulsePropagation.npz'
    oName="./FIGS/fig_optFrog_ESM_alpha%3.4lf.png"%(a0)

    def fetchData(fileLoc):
        data = np.load(fileLoc)
        return data['t'], data['Et']

    def windowFuncGauss(s0):
        return lambda t: np.exp(-t**2/2/s0/s0)/np.sqrt(2.*np.pi)/s0

    t,Et = fetchData(fName)

    res = optFrog(t,Et,windowFuncGauss,alpha=a0,tLim=(tMin,tMax,10), wLim=(wMin,wMax,3))
    spectrogramFigure((t,Et),res,oName=oName)


main()
# EOF: example_optFrog.py
