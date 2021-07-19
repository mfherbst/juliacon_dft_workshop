# Juliacon DFT workshop [![][binder-img]][binder-url]

These notebooks provide a brief mathematically-oriented introduction
into plane-wave density-functional theory (DFT)
and were prepared for the
[*A mathematical look at electronic structure theory*](https://pretalx.com/juliacon2021/talk/KK9KS7/)
workshop at Juliacon 2021.

The main topics of the notebooks are:
- Algorithms and numerical procedures to solve DFT (in particular SCF methods)
- Tools to understand and analyse convergence
- Opportunities to modern software engineering techniques
  (e.g. tracking floating-point error, automatic differentiation)
- Presentation of the [density-functional toolkit](https://dftk.org) (DFTK)
- Integration of DFTK within the Julia package ecosystem

For more details and to get started, see the
[Getting started](https://nbviewer.jupyter.org/github/mfherbst/juliacon_dft_workshop/blob/master/0_Getting_started.ipynb)
and [Installation](https://nbviewer.jupyter.org/github/mfherbst/juliacon_dft_workshop/blob/master/1_Installation.ipynb)
notebooks.

## Working with these notes online
If you do not want to install Julia, just run these notes
[on binder][binder-url],
which will allow you to play with the notebooks in your browser.
Note that for some of the exercises the computational performance
available on binder might not be sufficient.

[binder-url]: https://mybinder.org/v2/gh/mfherbst/juliacon_dft_workshop/master
[binder-img]: https://mybinder.org/badge_logo.svg
