#include <iomanip>
#include <iostream>
#include <string.h>
#include "Trans4dExample.h"
#include "Trans4d.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
//
using std::fixed;
using std::setprecision;

//initialize common i/o file stream objects
trans4d_example commonFiles;

/*  GetPoint
    Get a point from stdin via interaction from the user.
    User can specify whether to enter Lat/Long/Ellipsoid Height
    or cartesian X/Y/Z
    User also provides a name for this point up to 24 chars
*/
void GetPoint(int& LATD, int& LATM, double& SLAT, char& LATDIR,
              int& LOND, int& LONM, double &SLON, char& LONDIR, 
              char NAME24[24]/*string?*/,
              double& X, double& Y, double& Z,
              double& LAT, double& LON, double& EHT)
{
    DECLARE_COMMON_CONST
    double ELON;

    //retrieve the name of the point
    cout << " Enter name for point (24 character max)." << endl;
    string inputString;
    std::ws(cin); //skip whitespace characters
    std::getline(cin, inputString);
    strncpy(NAME24, inputString.c_str(), 24);
    
    label_select_coordinate_method:
    //retrieve the entry method for coordinates
    cout << " How do you wish to specify positional coordinates:" << endl;
    cout << "     1...geodetic latitude, longitude, ellipsoid height" << endl;
    cout << "     2...Cartesian (X,Y,Z) coordinates." << endl;

    int copt = 0;
    cin >> copt;
    if(copt == 1)
    {
        cout << " Enter latitude degrees-minutes-seconds in free format" << endl;
        cout << " with north being positive. For example,    35,17,28.3" << endl;
        cout << " For a point in the southern hemisphere, enter a minus sign" << endl;
        cout << " before each value. For example, -35,-17,-28.3" << endl;

        cin >> LATD >> LATM >> SLAT;

        cout << " Enter longitude degrees-minutes-seconds in free format" << endl;
        cout << " with west being positive.  To express a longitude measured" << endl;
        cout << " eastward, enter a minus sign before each value." << endl;

        cin >> LOND >> LONM >> SLON;
        
        cout << " Enter ellipsoid height in meters. (Note that" << endl;
        cout << " predicted motions are independent of this height.)" << endl;

        cin >> EHT;

        LAT =  (DBLE((LATD*60 + LATM)*60) + SLAT)/RHOSEC;
	    LATDIR = 'N';
        if(LAT < 0.e0)
        {
            LATD = -LATD;
            LATM = -LATM;
            SLAT = -SLAT;
            LATDIR = 'S';
        }
        LON = (DBLE((LOND*60 + LONM)*60) + SLON)/RHOSEC;
        ELON = -LON;
        trans4d::TOXYZ(LAT,ELON,EHT,X,Y,Z);
	    LONDIR = 'W';
	    if (LON < 0.0)
        {
            LOND = -LOND;
            LONM = -LONM;
            SLON = -SLON;
            LONDIR = 'E';
        }

    }
    else if (copt == 2)
    {
        cout << "Enter X coordinate in meters" << endl;
        cin >> X;

        cout << "Enter Y coordinate in meters" << endl;
        cin >> Y;

        cout << "Enter Z coordinate in meters" << endl;
        cin >> Z;
        if(!trans4d::FRMXYZ(X, Y, Z, LAT, LON, EHT))
            STOP(666);
        LON = -LON;
        if(LON < 0)
            LON = LON + TWOPI;
        int ISIGN;
        trans4d::TODMSS(LAT,LATD,LATM,SLAT,ISIGN);
	    LATDIR = 'N';
	    if (ISIGN == -1) 
            LATDIR = 'S';
        trans4d::TODMSS(LON,LOND,LONM,SLON,ISIGN);
	    LONDIR = 'W';
	    if (ISIGN == -1) 
            LONDIR = 'E';
    }
    else
    {
        cout << " Improper response -- try again." << endl;
        goto label_select_coordinate_method;
    }

    //todo C++ Port todo: error messages
    //       200 write (*,'(/)') 
    //       write (*,*) "Failed to read point name: ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   201 write (*,'(/)') 
    //       write (*,*) "Failed to read Coord. form option:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   202 write (*,'(/)') 
    //       write (*,*) "Failed to read latitude:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   203 write (*,'(/)') 
    //       write (*,*) "Failed to read longitude:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   204 write (*,'(/)') 
    //       write (*,*) "Failed to read ellipsoidal height:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   205 write (*,'(/)') 
    //       write (*,*) "Failed to read the X coordinate:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   206 write (*,'(/)') 
    //       write (*,*) "Failed to read the Y coordinate:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop

    //   207 write (*,'(/)') 
    //       write (*,*) "Failed to read the Z coordinate:ios=",ios
    //       write (*,*) "ABNORMAL TERMINATION"
    //       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
    //       stop
    
}

void Header()
{
    commonFiles.I2 << " Trans4D (VERSION " << TRANS4D_VERSION << ") OUTPUT" << endl;
    commonFiles.I2.flush(); //ensures any text still in the buffer is flushed output to file
}

//Interactively select a reference frame
void Menu1(int& kopt, string& mframe)
{
    int iframe[24+1];
    string nframe[24+1];

    iframe[1] = 1;
    nframe[1] = "NAD_83(2011/CORS96/2007)";   
    iframe[2] = 12;
    nframe[2] = "NAD_83(PA11/PACP00)";
    iframe[3] = 13;
    nframe[3] = "NAD_83(MA11/MARP00)";
    iframe[4] = 10;
    nframe[4] = "Stable NA (ITRF2014-PMM)";
    iframe[5] = 1;
    nframe[5] = "WGS_84(transit)";
//c     iframe[6] = 6 (This was incorrect in all versions of HTDP)
    iframe[6] = 5;
    nframe[6] = "WGS_84(G730)";
    iframe[7] = 8;
    nframe[7] = "WGS_84(G873)";
    iframe[8] = 11;
    nframe[8] = "WGS_84(G1150)";
    iframe[9] = 15;
    nframe[9]= "WGS_84(G1674)";
    iframe[10]= 15;
    nframe[10]= "WGS_84(G1762)";
    iframe[11]= 17;
    nframe[11]= "Pre-CATRF2022 =Caribbean";
    iframe[12]= 2;
    nframe[12]= "ITRF88";
    iframe[13]= 3;
    nframe[13]= "ITRF89";
    iframe[14]= 4;
    nframe[14]= "ITRF90/PNEOS_90/NEOS_90";
    iframe[15]= 5;
    nframe[15]= "ITRF91";
    iframe[16]= 6;
    nframe[16]= "ITRF92";
    iframe[17]= 7;
    nframe[17]= "ITRF93";
    iframe[18]= 8;
    nframe[18]= "ITRF94";
    iframe[19]= 8;
    nframe[19]= "ITRF96";
    iframe[20]= 9;
    nframe[20]= "ITRF97 or IGS97";
    iframe[21]= 11;
    nframe[21]= "ITRF2000 or IGS00/IGb00";
    iframe[22]= 14;
    nframe[22]= "ITRF2005 or IGS05";
    iframe[23]= 15;
    nframe[23]= "ITRF2008 or IGS08/IGb08";
    iframe[24]= 16;
    nframe[24]= "ITRF2014 or IGS14";

    cout << "  1...NAD_83(2011/CORS96/2007) (for use near North America) " << endl;
    cout << "  2...NAD_83(PA11/PACP00)      (for use on Pacific islands) " << endl;
    cout << "  3...NAD_83(MA11/MARP00)      (for use on the Mariana plate) " << endl;
    cout << "                                                   " << endl;
    cout << "  4...Relative to Stable North America according to         " << endl;
    cout << "      the ITRF2014 plate motion model        " << endl;
    cout << "                                                      " << endl;
    cout << "  5...WGS_84(transit) (NAD_83(2011) used)  15...ITRF91 " << endl;
    cout << "  6...WGS_84(G730) (ITRF91 used)           16...ITRF92 " << endl;
    cout << "  7...WGS_84(G873) (ITRF94 used)           17...ITRF93 " << endl;
    cout << "  8...WGS_84(G1150) (ITRF2000 used)        18...ITRF94 " << endl;
    cout << "  9...WGS_84(G1674) (ITRF2008 used)        19...ITRF96 " << endl;
    cout << " 10...WGS_84(G1762) (IGb08 used)   20...ITRF97 or IGS97" << endl;
    cout << " 11...Pre-CATRF2022 CARIBBEAN IFVM 21...ITRF2000 or IGS00/IGb00" << endl;
    cout << " 12...ITRF88                       22...ITRF2005 or IGS05 " << endl;
    cout << " 13...ITRF89                       23...ITRF2008 or IGS08/IGb08" << endl;
    cout << " 14...ITRF90 or (PNEOS_90/NEOS_90) 24...ITRF2014 or IGS14      " << endl;

    int iopt;
    std::ws(cin); //skip whitespace
    std::cin >> iopt;

    if(1 <= iopt && iopt <= 24)
    {
        mframe = nframe[iopt];
	    kopt = iframe[iopt];
    }
    else
    {
        mframe = mframe = "                ";
        kopt = iopt;
    }
}

void PrintProgramDescription()
{
    char buf[10];
    snprintf(buf, sizeof(buf), TRANS4D_VERSION);
    string version = buf;

    cout << " **************************************************" << endl;
    cout << " *  Trans4D (Transformations in 4 Dimensions)     *" << endl;
    cout << " *  SOFTWARE VERSION " << version                    << endl;
    cout << " *                                                *" << endl;
    cout << " *  AUTHORS:  R. Snay & C. Pearson & J. Saleh     *" << endl;
    cout << " *            Email: rssnay@aol.com               *" << endl;
    cout << " *                                                *" << endl;
    cout << " *  Fortran to C++ Port by C. Tewalt              *" << endl;
    cout << " **************************************************" << endl;
    cout << endl;
    cout << " This software incorporates numerical models that" << endl
         << " characterize continuous crustal motion as well as" << endl
         << " the episodic motion associated with earthquakes." << endl;

    cout << "The User Guide contains additional information and a set" << endl
         << " of exercises to familiarize users with the software." << endl;

    cout << " DISCLAIMER" << endl;
    cout << " The Trans4D software and supporting information are " << endl;
    cout << " currently distributed free of charge and are used by" << endl;
    cout << " the recipient with the understanding that the providers"<<endl;
    cout << " make no warranties, expressed or implied, concerning" << endl;
    cout << " the accuracy, completeness, reliabilty or suitability" << endl;
    cout << " of this software, of its constituent parts, or of any" << endl;
    cout << " supporting data." << endl;

    cout << " The providers shall be under no liability whatsoever" << endl;
    cout << " resulting from the use of this software. This software"<<endl;
    cout << " should not be relied upon as the sole basis for" << endl;
    cout << " solving a problem whose incorrect solution could" << endl;
    cout << " result in injury to person or property." << endl;
    cout << " Hit ENTER to continue." << endl;

    string line;
    std::getline(std::cin, line);
}

void ProgramLoop()
{
    while(true)
    {
        cout <<" ***************************************" << endl;
        cout <<" MAIN MENU:" << endl
            <<"    0... Exit software." << endl
            <<"    1... Estimate crustal velocities." << endl
            <<"    2... Estimate crustal displacements between dates." << endl
            <<"    3... Transform positions and/or observations," << endl
            <<"           entered in Blue Book format, across time" << endl
            <<"           and between reference frames." << endl
            <<"    4... Transform positions, entered in other formats," << endl
            <<"           across time and between reference frames." << endl
            <<"    5... Transform velocities between reference frames." << endl;

        int n;

        label_30_read_option:
        cin >> n;
        switch (n)
        {
        case 0:
            //exit program loop
            return;
        case  1:
            //estimate crustal velocities
            VELOC();
            break;
        case 2:
            // DPLACE();//todo
            break;
        case 3:
            cout << "  Unsupported option--select again" << endl;
            goto label_30_read_option;
            break;
        case 4:
            break;
        case 5:
            break;
        default:
            cout << "  Improper entry--select again" << endl;
            goto label_30_read_option;
            break;
        }
    }
}

void VELOC()
{
    int LATD, LATM, LOND, LONM;
    double SLAT, SLON, LAT, LON, EHT, X, Y, Z;
    char LATDIR, LONDIR;
    char NAME24[24];
    char PVOUT;
    string BLAB = "OUTSIDE OF REGION";

    double YLAT, YLON;
    double VN, VE, VU, VX, VY, VZ;
    double SN, SE, SU, SX, SY, SZ;
    int JREGN;

    cout << " Please enter name for the file to contain the" << endl
         << " predicted velocities." << endl;

    string NAMEF;
    std::ws(cin); //skip whitespace
    std::getline(cin, NAMEF);
    commonFiles.I2.open(NAMEF, std::fstream::in | std::fstream::out | std::fstream::app);

    Header();
    
    // Choosing reference system for velocities
    label_select_reference_frame:
    cout << "*********************************************" << endl;
    cout << " Select the reference frame to be used for specifying" << endl;
    cout << " positions and velocities." << endl;

    int iopt;
    string frame1;
    Menu1(iopt, frame1); //select reference frame

    if (iopt >= 1 && iopt <= numref)
    {
        commonFiles.I2 << "Velocities (with standard deviations) in mm/yr" << endl;
        commonFiles.I2 << "relative to " << frame1 << endl;
    }
    else
    {
        cout << "Improper selection -- try again" << endl;
        goto label_select_reference_frame;
    }
    
    // Choosing input format for locations where velocities are to be predicted
    cout << " ************************************************" << endl;
    cout << " Velocities will be predicted at each point whose" << endl;
    cout << " horizontal position is specified.  Please indicate"<< endl;
    cout << " how you wish to supply positions." << endl;

    cout << "    0... No more points. Return to main menu." << endl;
    cout << "    1... Individual points entered interactively." << endl;
    cout << "    2... Points on a specified grid." << endl;
    cout << "    3... The *80* records in a specified blue-book file." << endl;
    cout << "    4... Points on a specified line.  " << endl;
    cout << "    5... Batch file of delimited records of form: " << endl;
    cout << "         LAT,LON,TEXT " << endl;
    cout << "         LAT = latitude in degrees (positive north/DBL PREC) " << endl;
    cout << "         LON = longitude in degrees (positive west/DBL PREC) " << endl;
    cout << "         TEXT = Descriptive text (CHARACTER*24) " << endl;
    cout << "         Example:  " << endl;
    cout << "         40.731671553,112.212671753,SALT AIR " << endl;

    int OPTION;
    cin >> OPTION;

    if(OPTION == 0)
    {
        commonFiles.I2.close();
        if(PVOUT == 'Y') commonFiles.I3.close();
        return;
    }
    else if (OPTION == 1)
    {
        GetPoint(LATD, LATM, SLAT, LATDIR, LOND, LONM, SLON, LONDIR,
        NAME24, X, Y, Z, LAT, LON, EHT);
        trans4d::GTOVEL(LAT, LON, EHT, VN, VE, VU, VX, VY, VZ, JREGN,
            iopt, SN, SE, SU, SX, SY, SZ);
        
        if(JREGN == 0 )
        {
            cout << " ************************************* " << endl;
            cout << " A velocity can not be estimated because" << endl;
            cout << " the point is outside of the modeled region." << endl;
            cout << " For additional velocities, please indicate how" << endl;
            cout << " you wish to supply the horizontal coordinates." << endl;
        }
        else
        {  
// 	       WRITE(LUOUT,90) VN,SN,VE,SE,VU,SU,
//      &                        VX,SX,VY,SY,VZ,SZ
            cout << setprecision(2) << fixed;
            //todo: C++ port - nice column formatting
            cout << " **************************************" << endl;
            cout << " Northward velocity = " << VN << " +/- " << SN << " mm/yr" << endl;
            cout << " Eastward velocity  = " << VE << " +/- " << SE << " mm/yr" << endl;
            cout << " Upward velocity    = " << VU << " +/- " << SU << " mm/yr" << endl;

            cout << " X-dim. velocity    = " << VX << " +/- " << SX << " mm/yr" << endl;
            cout << " Y-dim. velocity    = " << VY << " +/- " << SY << " mm/yr" << endl;
            cout << " Z-dim velocity     = " << VZ << " +/- " << SZ << " mm/yr" << endl;

            cout << " **************************************" << endl;
            cout << " For additional velocities, please indicate how " << endl;
            cout << " you wish to specify the horizontal coordinates." << endl;

            //todo: C++ port write out to commonFiles.I2
// 	       WRITE(I2,1060)NAME24,LATD,LATM,SLAT,LATDIR,VN,SN,LOND,LONM, 
//      1                SLON, LONDIR,VE,SE,EHT,VU,SU
//                write(i2,1061) X,VX,SX,Y,VY,SY,Z,VZ,SZ
//  1060          FORMAT(/10X,A24,/
//      1'LATITUDE   = ',2I3,F9.5,1X,A1,'  NORTH VELOCITY =',F7.2,' +/- ',
//      & f4.2,/
//      2'LONGITUDE  = ',2I3,F9.5,1X,A1,'  EAST VELOCITY  =',F7.2,' +/- ',
//      & f4.2,/ 
//      3'ELLIPS. HT. = ',F10.3,' m', 6X,'UP VELOCITY    =',F7.2,' +/- ',
//      & f4.2)
//  1061  Format(
//      & 'X =',F13.3,' m',14X,'X VELOCITY     =',F7.2,' +/- ',
//      &    f4.2,/
//      5 'Y =',F13.3,' m',14X,'Y VELOCITY     =',F7.2,' +/- ',
//      &    f4.2,/
//      6 'Z =',F13.3,' m',14X,'Z VELOCITY     =',F7.2,' +/- ',
//      &    f4.2)
//   100          FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,F8.5,
//      1            1X,A1,1X,3(F8.2,' +/- ',F4.2))
//                IF(PVOUT .EQ. 'Y') THEN
//                    CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
//                ENDIF
//            ENDIF
        }
    }


}

int main()
{
    trans4d::InitBlockData();
    trans4d::MODEL();
    trans4d::SETTP();
    trans4d::SETRF();

    PrintProgramDescription();
    ProgramLoop();
    cout << endl;
}


