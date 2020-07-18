#include <iostream>
#include "trans4d.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;

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
    cin >> n;
    if(n == 0)
        return;
    
}

int main()
{
    trans4d trans4d;
    trans4d.InitBlockData();
    trans4d.MODEL();
    trans4d.SETTP();
    trans4d.SETRF();

    PrintProgramDescription();
    ProgramLoop();
    cout << endl;
}


