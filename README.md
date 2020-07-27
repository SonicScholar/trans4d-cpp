# TRANS4D C++

This is a port of TRANS4D (originally written in Fortran) into C++.
TRANS4D can be used to transform coordinates across reference frames and epochs.
It can also be used to predict velocities and displacements due to crustal mostion.

Note: TRANS4D is originally a Fortran program, so for ease of porting/translating,
many aspects of this port will feel like a conglomeration of Fortran and C/C++.

## Porting strategy
- Fortran arrays are indexed starting at 1, whereas C/C++ arrays are indexed starting at 0. 
To account for this, all arrays in this port start at 1 so as to minimize chances for
regression bugs.

- Use of goto/labels is supported in Fortran and C++, although the use of goto within C/C++
is generally discouraged as a best practice. However, for ease of porting, goto/labels
will still be used in a few places.

- Subroutine calls in Fortran are all "pass by reference" so the C++ translation from
subroutines to methods will involve passing variables by reference as well. 

- Fortran logical units for file i/o are replaced by std::fstream objects.

- Common block data are ported into structs with macros for declaring local reference
variables that reference the values in the structs. i.e. DECLARE_COMMON_XXXXXX

- Some char arrays will be replaced by std::string

- All scientific notation literal numeric values in C++ are interepreted as 8 byte doubles.

- UI related code is bein separated into the main sample exe while core functionality is
refactored into a shared library (dll in Windows)

- C++ has no implicit variable types. In cases where some variables were arbitrarily prefixed with
X or Y to ensure they were doubles, those prefixes have been removed to avoid confusion. For example,
the variable LON might be desired as a double, but implicitly be an int due to the following declaration:

    ```
    IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    IMPLICIT INTEGER*4 (I-N)
    ```

- Bluebook format related operations are unsupported at this time.

# Compile
```
mkdir build
meson --prefix=<PATH>/bin build
ninja -C build install
```

or...

```
meson builddir
cd builddir
meson compile
```

# Current Porting Progress

All of the subroutines are listed here from the original trans4d.f file.
Below shows progress on subroutines being ported, as well as what is complete.
The `In Progress` section shows what line number the subroutine is translated up to.

|Key| Status     |
|---|:-----------|
| X | Not needed |
| * | In progress|
| ? | Not started|


## IN PROGRESS

| SUBROUTINE     | STATUS         | Progress       |
|----------------|:--------------:|----------------|
| TRANS4D MAIN   | *              | line 124. skipped line 122: call to dplace |
| VELOC          | *              | line 2284. Only option1 at this time. skipped writing to file      |





## COMPLETE

| SUBROUTINE              | STATUS         | Notes          |
|-------------------------|:--------------:|----------------|
| COMVEL                  | Done           |
| FRMXYZ                  | Done           |
| GETBDY                  | Done           |
| GETGRID                 | Done           |
| GETPNT                  | Done           |
| GETREG                  | Done           |
| GRDVEC                  | Done           |
| GRDWEI                  | Done           |
| GTOVEL                  | Done           |
| HEADER                  | Done           |
| IYMDMJ                  | Done           |
| MENU1                   | Done           |
| MODEL                   | Done           |
| PLATVL                  | Done           |
| POLYIN                  | Done           | 
| RADII                   | Done           |
| SETRF                   | Done           |
| SETTP                   | Done           |
| to_itrf2014             | Done           |
| TODMSS                  | Done           | 
| TOXYZ                   | Done           |
| TOVNEU                  | Done           |
| TOVXYZ                  | Done           |
| VTRANF                  | Done           |
| XTOITRF2014             | Done           |
| to_std_dev_xyz_velocity | Done           |


## TODO

| SUBROUTINE     | STATUS         | NOTES          |
|----------------|:--------------:|----------------|
| CHECK          | ?              |
| COMPSN         | ?              |
| DDXYZ          | ?              |
| DIRCT1         | ?              |
| DISLOC         | ?              |
| DPLACE         | ?              |
| DSDA           | ?              |
| extract_name   | ?              |
| from_itrf2014  | ?              |
| get_frames     | ?              |
| GETLYN         | ?              |
| GETMDY         | ?              |
| GETPO4         | ?              |
| GETVLY         | ?              |
| GRDAMP         | ?              |
| GRDCHK         | ?              |
| HELINV         | ?              |
| interprate_latlon_record        | ?              |  refactor spelling?
| interprate_latlonvel_record     | ?              |  refactor spelling?
| interprate_XYZ_record           | ?              |  refactor spelling?
| NEWCOR         | ?              |
| OKADA          | ?              |
| OKADAW         | ?              |
| PDISP          | ?              |
| PREDV          | ?              |
| PRNTTP         | ?              |
| PRINTVL        | ?              |
| PSGWEI         | ?              |
| PVPRNT         | ?              |
| RADR8T         | ?              |
| RDEG           | ?              |
| RFCON          | ?              |
| RFCON1         | ?              |
| TOCHAR         | ?              |
| TOMNT          | ?              |
| TNFDAT         | ?              |
| tran_frames    | ?              |
| TRAVEC         | ?              |
| trfbb          | ?              |
| TRFDAT         | ?              |
| TRFPOS1        | ?              |
| TRFVEL         | ?              |
| UPBB4          | ?              |
| UPGFI4         | ?              |

---------------
