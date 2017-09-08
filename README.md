# Sulley

THIS PROJECT IS NO LONGER ACTIVE. GO TO https://gitlab.cern.ch/alvarezc/Sulley

Fitting monster.

The main classes are:

* toystudy1DHistBremCatExpoBkg performs a toyMC study using a 1D fit
* toystudy2DHistBremCatTsallisBkg performs a toyMC study using a 2D fit
* toystudy performs a toyMC study (1D/2D/constrained)

Tested with ROOT 6.04/02

The HistFactory implementation of the fit needs a patched version of ROOT that fixes the implementation of multidimensional fits. This patched is installed in /vols/build/lhcb/ROOT_6.06.02_patch_histfactory at the lx0y machines. To return to the default version of ROOT, comment out the relevant lines in setup.sh.

INSTALL

Setup the environment for ROOT

$ source setup.sh

$ make

