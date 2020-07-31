# TRANS4D C++

This is a port of TRANS4D (originally written in Fortran) into C++.
TRANS4D can be used to transform coordinates across reference frames and epochs.
It can also be used to predict velocities and displacements due to crustal mostion.

## About

This is currently a port of version 0.2.6 of TRANS4D. The C++ port is a volunteer project
made freely available to help software engineers to be able to compile HTDP/Trans4d on
more devices than can be targeted with Fortran. Before using this for production purposes,
please make sure to test coordinates using the original TRANS4D program. At the time of
writing, many test cases need to be written to verify the accuracy of the translation of
velocity grid, earthquake, and post seismic activity.



## TRANS4D vs HTDP
The National Geodetic Survey's (NGS) official tool for calculating position transformations
and crustal motion is Horizontal Time Dependent Positioning (HTDP). TRANS4D is the evolution
of that program, written and maintained by Dr. Richard Snay (NGS, Retired). He continues
to improve the resolution of velocity grids and incorporates new models into his software.

As of mid 2020, Trans4d's models are more precise than HTDP with velocity grids containing 
points with 6 km resolution in CONUS. Trans4d's models also incorporate vertical crustal
motion. In conversations with Dr. Snay and several NGS employees maintaining HTDP, the NGS
plans to use Trans4d to upgrade the functionality of HTDP.

## Transforming Positions

- Make sure the Data4.2.5A.txt and Data4.2.5B.txt files are in the same directory as the
executable

- Include "Trans4d.h"

- Use the `trans4d::TransformPosition() method`

- See the table below for values to use in the IOPT parameter
    ```// 1;  "NAD_83(2011/CORS96/2007)";   
    // 1;  "WGS_84(transit)";
    // 10; "Stable NA (ITRF2014-PMM)";
    // 12; "NAD_83(PA11/PACP00)";
    // 13; "NAD_83(MA11/MARP00)";
    // /***iframe[6] = 6 (This was incorrect in all versions of HTDP)***/
    // 5;  "WGS_84(G730)";
    // 8;  "WGS_84(G873)";
    // 11; "WGS_84(G1150)";
    // 15; "WGS_84(G1674)";
    // 15; "WGS_84(G1762)";
    // 17; "Pre-CATRF2022 =Caribbean";
    // 2;  "ITRF88";
    // 3;  "ITRF89";
    // 4;  "ITRF90/PNEOS_90/NEOS_90";
    // 5;  "ITRF91";
    // 6;  "ITRF92";
    // 7;  "ITRF93";
    // 8;  "ITRF94";
    // 8;  "ITRF96";
    // 9;  "ITRF97 or IGS97";
    // 11; "ITRF2000 or IGS00/IGb00";
    // 14; "ITRF2005 or IGS05";
    // 15; "ITRF2008 or IGS08/IGb08";
    // 16; "ITRF2014 or IGS14";
    ```
    Example code converting from NAD83(2011) epoch 2010.0 to ITRF2014 epoch 2020.0
    
    ```double latDegrees = 40.0001;
    double lonDegrees = -105.0001;
    double eht = 1500;

    double inDate = 2010.0;
    double outDate = 2020.0;

    //todo: make these values constants or enums
    int nad83Opt = 1;
    int itrf2014Opt = 16;

    double newLat, newLon, newEht;
    trans4d::TransformPosition(latDegrees, lonDegrees, eht, nad83Opt, itrf2014Opt, inDate, outDate, newLat, newLon, newEht);
    ```

## Porting strategy
TRANS4D is originally a Fortran program, so for ease of porting/translating,
many aspects of this port will feel like a conglomeration of Fortran and C/C++.

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

## Contact

For any inquiries on this software:
- feel free to reach out to me at `collin@sonicscholar.com`
- Dr. Snay's email is `rssnay@aol.com`
- For inquiries about HTDP, contact the NGS at `ngs.cors.htdp@noaa.gov`
