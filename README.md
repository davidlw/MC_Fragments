To be run in the examples directory of PYTHIA 8.310 by replacing the existing makefile

It is expected you installed pythia with the ROOT and fastjet additional packages

The 'new' configurations are:

pp_highMultGen -> LHC jets with potential to filter on jet multiplicity
DIS -> lepton-proton (for EIC or DIS physics)
LEP1 -> lepton-lepton at Z pole (only Z->qqbar process for now)
photonHadron -> anything with a photon hitting a proton (photoproduction or UPC - see code for more details)

You can build the relevant config with, for example:

make DIS

and run with

./DIS <jobnumber>


