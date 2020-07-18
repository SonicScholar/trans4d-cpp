# TRANS4D C++

This is a port of TRANS4D (originally written in Fortran) into C++.
TRANS4D can be used to transform coordinates across reference frames and epochs.
It can also be used to predict velocities and displacements due to crustal mostion.

Note: TRANS4D is originally a Fortran program, so for ease of porting/translating,
many aspects of this port will feel like a conglomeration of Fortran and C/C++.

-----------------------
## Current Porting Progress

_Key_

- X - Not needed
- \* - In progress
- ? - Not started
- D - Done

### IN PROGRESS

* PROGRAM TRANS4D
    - line 101: todo: figure out if(ios /= 0) goto 51
	- line 117:



### TODO

- ? GETREG
- ? POLYIN
- ? RADR8T
- ? COMVEL
- ? getgrid
- ? PLATVL
- ? PVPRNT
- ? DSDA
- ? HELINV
- ? TRFDAT
- ? TNFDAT
- ? TOXYZ
- ? RADII
- ? GETLYN
- ? DIRCT1
- ? TOCHAR
- ? DDXYZ
- ? DISLOC
- ? OKADA
- ? OKADAW
- ? GRDWEI
- ? GRDVEC
- ? RDEG
- ? TOMNT
- ? TOVNEU
- ?  TOVXYZ
- ? to_std_dev_xyz_velocity
- ? DPLACE
- ? VELOC
- ? GTOVEL
- ? XTOITRF2014
- ? TRFPOS1
- ? PRNTTP
- ? from_itrf2014
- ? to_itrf2014
- ? menu1
- ? trfbb
- ? GETPO4
- ? UPBB4
- ? CHECK
- ? UPGFI4
- ? TODMSS
- ? RFCON
- ? RFCON1
- ? TRAVEC
- ? GETPNT
- ? GETVLY
- ? COMPSN
- ?  NEWCOR
- ? PREDV
- ? TRFVEL
- ? VTRANF
- ? PRINTVL
- ? HEADER
- ? GETMDY
- ? PDISP
- ? GRDCHK
- ? PSGWEI
- ? GRDAMP
- ? extract_name
- ? interprate_XYZ_record       - refactor spelling?
- ? interprate_latlon_record    - refactor spelling?
- ? interprate_latlonvel_record - refactor spelling?
- ? get_frames
- ? tran_frames

### COMPLETE
--------------
- D MODEL
- D GETBDY
- D SETTP
- D SETRF
- ? IYMDMJ