CoNSEnsX
========

You find the current working (Perl) version [here](http://himalia.chem.elte.hu/cgi-bin/consensx.cgi).

## New functions:

#### I/O modifications:

##### Client side
* Completely rewritten HTML5 client side
* JavaScript form checking in client side
* PDB file is automaticly downdloaded if PDB ID specified on client side
* Optional PDB model alignment
    * All models will be fitted to the first model with an RMSD minimazing transformation
    * Fitting range is now selecteble on client side

##### CLI
* BMBR restraint file format for NOE distances

##### HTML output
* All calculation data goes to a working directory, each calculation has its own ID and HTML page
* All images are now scalable (SVG)

## Implemented back-calculations:

##### Structural parameters
* NOE distance restraint violations
    * NOE distances are calculated and averaged for the whole ensemble
    * Count of total violations
    * The distribution of violation distances is visualized by histogram
    * [new] Optional r-6 averaging

* PRIDE-NMR
    * Scores similarity between the NOE pattern and the back-calculated H-H distance distribution
    * Model with best/worst score
    * Average score of the models and standard deviation
    * Informative diagram


##### Dynamical parameters
* Residual dipolar couplings (RDC)
    * Correspondence of experimental RDC data with the alignments predicted by PALES (Prediction of ALignmEnt from Structure)
    * RDC data sets are now aligned separately (one list -> one alignment)
    * separate alignments have separate section in HTML output

* Scalar coupling
    * User can select the used set of Karplus equation constant

* Chemical shifts
    * Correspondence of experimental chemical shifts data with the chemical shifts predicted by SHIFTX

* Order parameter (S<sup>2</sup>)
    * Correspondence of experimental order parameter data with order parameters back-calculated from the input PDB models
    * Only N-H order parameter is implemented
