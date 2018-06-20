# optfrog 

> Analytic signal spectrograms with optimized time-frequency resolution

The Python package optfrog comprises tools for the calculation of spectrograms
with optimized time and frequency resolution. Gabor's uncertainty
principle prevents both from beeing optimal simultaneously for a given window
function employed in the underlying short-time Fourier analysis.  
Aim of optfrog is to provide a time-frequency representation of the input
signal with marginals that represent the original intensities per unit time and
frequency similarly well. A tunable paramter might be adjusted to emphasize
time or frequency resolution.  


### Prerequisites

Optfrog requires the functionality of 

* Numpy
* Scipy
* Matplotlib


### Installing

Get a local copy of the repository and install optfrog by running

```
python setup.py install
```

from the commandline.

In case you downloaded the source distribution optfrog-1.0.0.tar.gz from the folder dist onto a Unix system with user rights only, run

```
tar -xvf 
cd optfrog-1.0.0
python setup.py install --user
```

from the commandline. This, however will only install the optfrog tools without the examples folder.

## Example Programs

Optfrog also comes with several sample programs in the examples directory. For example, 
changing to the examples folder

```
python example_optFrog.py
```

will compute a time-frequency optimized spectrogram for an input signal characterizing 
a short intense optical pulse after propagation in a nonlinear waveguide.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

This work received funding from the VolkswagenStiftung within the
‘Niedersächsisches Vorab’ program in the framework of the project ‘Hybrid
Numerical Optics’ (HYMNOS; Grant ZN 3061). 

