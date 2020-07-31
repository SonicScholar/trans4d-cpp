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






## COMPLETE

| SUBROUTINE              | STATUS         | Notes          |
|-------------------------|:--------------:|----------------|
| COMVEL                  | Done           |
| COMPSN                  | Done           |
| DISLOC                  | Done           |
| FRMXYZ                  | Done           |
| from_itrf2014           | Done           |
| GETBDY                  | Done           |
| GETGRID                 | Done           |
| GETREG                  | Done           |
| GETPNT                  | Done           |
| GRDAMP                  | Done           |
| GRDCHK                  | Done           |
| GRDVEC                  | Done           |
| GRDWEI                  | Done           |
| GTOVEL                  | Done           |
| HEADER                  | Done           |
| IPSGRD                  | Done           |
| IUNGRD                  | Done           |
| IYMDMJ                  | Done           |
| MENU1                   | Done           |
| MODEL                   | Done           |
| NEWCOR                  | Done           |
| OKADA                   | Done           |
| OKADAW                  | Done           |
| PLATVL                  | Done           |
| POLYIN                  | Done           | 
| PREDV                   | Done           |
| PSDISP                  | Done           |
| PSGWEI                  | Done           |
| RADII                   | Done           |
| RADR8T                  | Done           |
| SETRF                   | Done           |
| SETTP                   | Done           |
| to_itrf2014             | Done           |
| to_std_dev_xyz_velocity | Done           |
| TODMSS                  | Done           | 
| TOMNT                   | Done           |
| TOVNEU                  | Done           |
| TOVXYZ                  | Done           |
| TransformPosition       | Done           | Adaptation of TRFPOS1
| TOXYZ                   | Done           |
| VELOC                   | Done           | Only option1 at this time. |
| VTRANF                  | Done           |
| XTOITRF2014             | Done           |



## TODO

| SUBROUTINE     | STATUS         | NOTES          |
|----------------|:--------------:|----------------|
| CHECK          | ?              |
| DDXYZ          | ?              |
| DIRCT1         | ?              | Nice Utility function
| DPLACE         | ?              |
| DSDA           | ?              |
| extract_name   | X              |
| get_frames     | X              |
| GETLYN         | X              |
| GETMDY         | X              |
| GETPO4         | X              |
| GETVLY         | X              |
| HELINV         | ?              | Nice Utility function
| interprate_latlon_record        | ?              |  refactor spelling?
| interprate_latlonvel_record     | ?              |  refactor spelling?
| interprate_XYZ_record           | ?              |  refactor spelling?
| PRNTTP         | X              |
| PRNTVL         | X              |
| PVPRNT         | X              |
| RDEG           | ?              | Nice Utility function
| RFCON          | ?              |
| RFCON1         | ?              |
| TOCHAR         | X              |
| TNFDAT         | X              |
| tran_frames    | ?              | Not originally part of HTDP. No crustal motion accounted for here.
| TRAVEC         | ?              | Nice Utility function
| trfbb          | X              |
| TRFDAT         | X              |
| TRFPOS1        | X              | Captured this functionality in TransformPosition
| TRFVEL         | ?              | Nice Utility Function
| UPBB4          | X              |
| UPGFI4         | X              |

---------------