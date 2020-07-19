# TRANS4D C++

This is a port of TRANS4D (originally written in Fortran) into C++.
TRANS4D can be used to transform coordinates across reference frames and epochs.
It can also be used to predict velocities and displacements due to crustal mostion.

Note: TRANS4D is originally a Fortran program, so for ease of porting/translating,
many aspects of this port will feel like a conglomeration of Fortran and C/C++.

## Porting strategy
Fortran arrays are indexed starting at 1, whereas C/C++ arrays are indexed starting at 0. 
To account for this, all arrays in this port start at 1 so as to minimize chances for
regression bugs.

Use of goto/labels is supported in Fortran and C++, although the use of goto within C/C++
is generally discouraged as a best practice. However, for ease of porting, goto/labels
will still be used in a few places.

Subroutine calls in Fortran are all "pass by reference" so the C++ translation from
subroutines to methods will involve passing variables by reference as well. 

Fortran logical units for file i/o are replaced by std::fstream objects.

Common block data are ported into structs with macros for declaring local reference
variables that reference the values in the structs. i.e. DECLARE_COMMON_XXXXXX

Some char arrays will be replaced by std::string

All scientific notation literal numeric values in C++ are interepreted as 8 byte doubles.

UI related code is bein separated into the main sample exe while core functionality is
refactored into a shared library (dll in Windows)

Bluebook format related operations are unsupported at this time.


-----------------------
## Current Porting Progress

All of the subroutines are listed here from the original trans4d.f file.
Below shows progress on subroutines being ported, as well as what is complete.
The `In Progress` section shows what line number the subroutine is translated to.


_Key_

- X - Not needed
- \* - In progress
- ? - Not started
- D - Done

### IN PROGRESS
--------

* PROGRAM TRANS4D
	- line 124: estimate crustal velocities
	- (skipped line 122: call to dplace)

* VELOC
	- 2236

* GETPNT
	- 4756



### TODO
--------
? DPLACE
? GETREG
? POLYIN
? RADR8T
? COMVEL
? getgrid
? PLATVL
? PVPRNT
? DSDA
? HELINV
? TRFDAT
? TNFDAT
? TOXYZ
? RADII
? GETLYN
? DIRCT1
? TOCHAR
? DDXYZ
? DISLOC
? OKADA
? OKADAW
? GRDWEI
? GRDVEC
? RDEG
? TOMNT
? TOVNEU
? TOVXYZ
? to_std_dev_xyz_velocity
? VELOC
? GTOVEL
? XTOITRF2014
? TRFPOS1
? PRNTTP
? from_itrf2014
? to_itrf2014
? trfbb
? GETPO4
? UPBB4
? CHECK
? UPGFI4
? TODMSS
? RFCON
? RFCON1
? TRAVEC
? GETPNT
? GETVLY
? COMPSN
? NEWCOR
? PREDV
? TRFVEL
? VTRANF
? PRINTVL
? GETMDY
? IYMDMJ
? PDISP
? GRDCHK
? PSGWEI
? GRDAMP
? extract_name
? interprate_XYZ_record       - refactor spelling?
? interprate_latlon_record    - refactor spelling?
? interprate_latlonvel_record - refactor spelling?
? get_frames
? tran_frames

### COMPLETE
--------------
D MENU1
D MODEL
D GETBDY
D HEADER
D SETTP
D SETRF