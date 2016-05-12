#ifndef _MC_H
#define _MC_H

#include "alloc.h"
#include "Stream.h"

/************************************************************
*** MC simulation using various ensembles:
    NVT
    NPT
    NPH (t = 0)
    NtH
*** References
    [For NtH]
    J. R. Ray et al, J. Chem. Phys., 1984
    P. J. Fay et al, PRA, 1992
    M. Karimi et al, PRB, 1998
*** Parameters:
    seed         --> seed for SFMT RNG
    mol          --> initial molecule
    enthalpy     --> fixed enthalpy [kT = 2(H-U(r))/(3N)]
    dxmax        --> maximum displacement for an atom's x/y/z
    dhmax        --> maximum displacement of the cell's a/b/c
    rStep        --> MC steps in relaxation (1 step = N trial moves)
    eStep        --> MC steps in equilibrium (recorded statistics)
    hStep        --> try moving h(cell) every "hStep" MC steps
    stream       --> deal with file io
    NOTICE: file dump only happens in the equilibrium stage
*************************************************************/
namespace MC {

void useNVT(long seed,
	    AMod::Molecule& mol,
	    MCBasic& mc,
	    MCAgent& mcAgent,
	    double kT,
	    double dxmax,
	    long rStep,
	    long eStep,
	    Stream& stream);

void useNPT(long seed,
	    AMod::Molecule& mol,
	    MCBasic& mc,
	    MCAgent& mcAgent,
	    double pressure,
	    double kT,
	    double dxmax,
	    double dhmax,
	    long rStep,
	    long eStep,
	    long hStep,
	    Stream& stream);

void useNPH(long seed,
	    AMod::Molecule& mol,
	    MCBasic& mc,
	    MCAgent& mcAgent,
	    double pressure,
	    double enthalpy,
	    double dxmax,
	    double dhmax,
	    long rStep,
	    long eStep,
	    long hStep,
	    Stream& stream);

void useNtH(long seed,
	    AMod::Molecule& mol,
	    MCBasic& mc,
	    MCAgent& mcAgent,
	    double pressure,
	    const double tension[9], //tension matrix (column major)
	    double enthalpy,
	    double volume0,          //volume0 = det(h0)
	    const double h0[9],      //h0 matrix (column major)
	    double dxmax,
	    double dhmax,
	    long rStep,
	    long eStep,
	    long hStep,
	    Stream& stream);

}/* MC */

#endif
