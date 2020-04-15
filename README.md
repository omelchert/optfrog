# optfrog 

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

The Python package optfrog comprises tools for the calculation of spectrograms
with optimized time and frequency resolution. Gabor's uncertainty
principle prevents both from beeing optimal simultaneously for a given window
function employed in the underlying short-time Fourier analysis.  
Aim of optfrog is to provide a time-frequency representation of the input
signal with marginals that represent the original intensities per unit time and
frequency similarly well. A tunable paramter might be adjusted to emphasize
time or frequency resolution.  


### Prerequisites

The tools provided by the `optfrog` package require the functionality of 

* numpy (>=1.8.0rc1)
* scipy (>=0.13.0b1)

Further, the sample scripts provided in the `examples` folder require the funcality of

* matplotlib (>=1.2.1)

for generating figures as the one shown below.

### Installing

Get a local copy of the repository and install optfrog by running

```
python setup.py install
```

from the commandline. The installation requires pythons setuptools.

In case you downloaded the source distribution `optfrog-1.0.0.tar.gz` from the folder `dist` onto a Unix system with user rights only, run

```
tar -xvf optfrog-1.0.0.tar.gz
cd optfrog-1.0.0/
python setup.py install --user
```

from the commandline. This, however will only install the optfrog tools without the examples folder. So make sure to also fetch a local copy of the examples provided here.

## Example Programs

`optfrog` also comes with several sample programs in the examples directory. For example, 
changing to the `examples` folder and running

```
python example_optFrog.py
```

will compute a time-frequency resolution optimized spectrogram for an input signal characterizing 
a short intense optical pulse after propagation in a nonlinear waveguide. The above example script will generate the below figure in subfolder `FIGS`.

![alt text](https://github.com/omelchert/optfrog/blob/master/examples/FIGS/fig_optFrog_ESM_alpha0.0000.png)

### Brief explanation of the above figure

The center part of the figure shows an `optFrog` spectrogram for an exemplary input signal. The intesity of the
spectrogram is normalized to a maximum value of unity. The subfigure on top allows to compare the intensity per unit time of the normalized analytic signal for the input data (gray line) to the time marginal obtained from the spectrogram (black line). The subfigure on the right shows the intensity per unit frequency of the analytic signal (gray line) as well as the frequency marginal computed from the spectrogram (black line).

## Availability of the software

We further prepared an [optFROG compute capsule](https://codeocean.com/capsule/7823161/tree) on [Code Ocean](https://codeocean.com), allowing to directly run and modify an exemplary simulation without the need to locally install the software. 

## Links

The optfrog software package is derived from our research software and is described in 

> O. Melchert, U. Morgner, B. Roth, and A. Demircan, "OptFROG — Analytic signal spectrograms with optimized time–frequency resolution", [SoftwareX 10, 100275 (2019)](https://doi.org/10.1016/j.softx.2019.100275)

Further use-cases demonstrating parts of the functionality of `optfrog` as a tool for signal analysis in ultrashort pulse propagation can be found under

> O. Melchert, U. Morgner, B. Roth, I. Babushkin, and A. Demircan, "Accurate propagation of ultrashort pulses in nonlinear waveguides using propagation models for the analytic signal," [Proc. SPIE 10694, Computational Optics II, 106940M (2018)](https://doi.org/10.1117/12.2313255)

and

> O. Melchert, S. Willms, S. Bose, A. Yulin, B. Roth, F. Mitschke, U. Morgner, I. Babushkin, and A. Demircan, "Soliton Molecules with Two Frequencies", [Phys. Rev. Lett. 123, 243905 (2019)](https://doi.org/10.1103/PhysRevLett.123.243905)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

This work received funding from the Deutsche Forschungsgemeinschaft  (DFG) under
Germany’s Excellence Strategy within the Cluster of Excellence PhoenixD
(Photonics, Optics, and Engineering – Innovation Across Disciplines) (EXC 2122,
projectID 390833453).
