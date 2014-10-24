CoNSEnsX
========

CoNSEnsX Python rewrite

You find the current working (Perl) version [here](http://himalia.chem.elte.hu/cgi-bin/consensx.cgi).

## Implemented functions:

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

#### Implemented back-calculations:

* RDC
    * RDC data sets are now aligned separately (one list -> one alignment)
    * separate alignments have separate section in HTML output
* Scalar coupling
* Chemical shifts

* Order parameter (S<sup>2</sup>)
    * Only N-H order parameter is implemented
