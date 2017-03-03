# MUltipactor in Lower Hybrid antenna waveguides (MULH)

This code does a power sweep to calculate at which point the multipactor
starts to develop inside a single rectangular waveguide.
EM fields can be resolved analytically for the TE10 mode of a rectangular
waveguide, or through the FDTD method with a TE10 excited at one end of the guide, or other
fields could be imported (e.g. from COMSOL).

The particle wall interaction is modeled using G. Chen, L. Liu (2011)
ejection-collection algorithm but the actual method for deciding between
absorption, elastic or inelastic collision or secondaries has been slightly
modified to match FEST3D results. The option to use the Furman & Pivi model
is also made available.

The Ecuyer Taus random number generator from Alan Miller's library is used.
It has a period of 2^88 and the poisson deviates are generated using this
generator and the waiting time method for means < 15.

Studies of how sensitive the threshold is to different parameters can be
performed through the control of the variable atype. At the moment only 11
different studies can be performed but it is not hard to program more into
MULH.

NOTE: There are some outstanding issues that ought to be resolved:
- MULH predicts a breakdown for wide waveguides (WR-229 and WR-284), while FEST3D doesn't.
The breakdown calculated by both programs for other geometries seem to agree though. Perhaps
the fields and particle steppers in these wide guides should be checked.
- MULH didn't see a resonance when applaying a DC toroidal B-field equal to
the one predicted by a gyrofrequency of 3.7GHz (Tore-Supra), while Spark3D did.
- The trends in breakdown when DC B-fields are applied were never rigorously checked against
FEST3D or Spark3D.

The MULH parent directory contains four subfolders.
 - sources: contains all the source files used in the code. subdirectory
	active contains the latest in-use source files.
 - data: .txt files containing data used in some runs and the subdirectory
	results, where output data is stored. It also contains two m-files
	used for processing COMSOL fields from Melanie Preyna's PICCOLO and
	for processing MULH output data.
 - execs: where the executables are placed during compilation
 - test: for storing executables and other files used while testing pieces
	of or all of the code.

To run the code go to the directory where the makefile is (execs) in the
terminal.
Then run the following commands:
```
make
export LD_LIBRARY_PATH="path to executable"
./MULH
```

The code was developed and worked well with the gfortran compiler.
It hasn't been tried with other compilers.

For more details on the code see:
Francisquez, M. 2012. Power limit modeling of lower hybrid antenna waveguides
in tokamaks. BA Honors Thesis. Dartmouth College, November 2012.

Please feel free to contact me with questions and/or suggestions.

AUTHOR: Manaure Francisquez
Department of Physics & Astronomy
Dartmouth College
Hanover, NH 03755
USA
http://engineering.dartmouth.edu/~d24789f

Spring 2011 - Spring 2012
