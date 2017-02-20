#ifndef __InitTree__
#define __InitTree__

#define GENE_DEBUG

#include "TTree.h"

using namespace std;

// event information =====================================================

extern int   eRun;     // run number
extern int   eSpill;   // spill number
extern int   eEvinspill;   // event number in spill
extern int   eTrig;    // trigger mask

extern double  eImu;   // incident muon intensity.
extern double  ePolu;   // polar upstream
extern double  ePold;   // polar downstream 
extern double  eD;      // depolarization factor

// vertex ================================================================
extern double vX;    // X coordinate of vertex
extern double vY;    // Y coordinate of vertex
extern double vZ;    // Z coordinate of vertex
extern int    vN;    // Number of tracks in vertex
extern double vChi2; // Chi2 of primary vertex

// MC block
extern double vmcX;    // X coordinate of MC vertex
extern double vmcY;    // Y coordinate of MC vertex
extern double vmcZ;    // Z coordinate of MC vertex
extern int    vmcN;    // Number of tracks in MC vertex

// kinematics ============================================================
extern double kEgamma; // energy of virtual photon
extern double kQ2;     //  - mass of virtual photon
extern double kXbj;    // bjorken variable
extern double kY;      // mu relative energy loss  
extern double kW2;     // W^2 varlable

// MC block
extern double kmcEgamma; // energy of virtual photon MC
extern double kmcQ2;     //  - mass of virtual photon MC
extern double kmcXbj;    // bjorken variable MC
extern double kmcY;      // mu relative energy loss MC  
extern double kmcW2;     // W^2 varlable MC

// mu,mup ================================================================

extern double bP;  // Momentum of beam
extern double bPx; // X mom component
extern double bPy; // Y mom component
extern double bE;  // Energy of beam
extern double mP;  // Momentum scattered muon
extern double mPx; // X mom component
extern double mPy; // Y mom component
extern double mE;  // Energy of mu'

// MC block
extern double bmcP; // Momentum of beam MC
extern double mmcP; // Momentum of mu' MC

extern double mZlast; // Z of mu' last hit
extern double mNhits; // # of hits in mu' track
extern double mTime;  // time of mu' track
extern double mX0;    // radiative length traversed by mu'


// hadrons ===============================================================
  
extern int hN;         // number of hadrons
extern double hMinv;   // invariant mass of hadron pair
extern double *hChi2;  // chi^2 of a hadron track
extern double *hQ;     // charge of hadron
extern double *hX_MM1V;
extern double *hY_MM1V;
extern double *hX_MM1U;
extern double *hY_MM1U;
extern double *hX_MM1Y;
extern double *hY_MM1Y;
extern double *hX_DC1;
extern double *hY_DC1;


// lab 
extern double *hP;     // hadron momentum
extern double *hE;     // hadron energy
extern double *hPx;    // X mom component
extern double *hPy;    // Y mom component
extern double *hPz;    // Z mom component
extern double *hPt;    // transverse momentum of hadron
extern double *hPl;    // longitudinal momentum of hadron
extern double *hZ;     // fraction of gamma energy carried by hadron
extern double *hPhi;   // azimuthal angle of hadron track
extern double *hTheta; // polar angle of hadron track

// cms
extern double *hXf;    // x Feynman variable of a hadron
  
// ID
extern double *hZlast; // Z of last hit in hadron track
extern double *hNhits; // # of hits in hadron track
extern double *hTime;  // time of hadron track
extern double *hX0;    // radiative lenght traversed by hadron

// calorimeter information
extern double *hEc;    // Energy deposition in all calos
extern double *hHc1;   // HC01
extern double *hHc2;   // HC02
extern double *hEc1;   // EC01
//L.S.
extern double *hEc2;   // EC02
extern double *hXCClus;
extern double *hYCClus;
extern double *hZCClus;


// MC block 

extern double *hmcP;     // hadron momentum MC
extern double *hmcPt;    // transverse momentum of hadron MC
extern double *hmcTheta; // polar angle of hadron track MC
extern double *hmcZ;     // fraction of gamma energy carried by hadron MC
extern double *hmcXf;    // x Feynman variable of a hadron MC
extern int    *hmcPid;   // Particle ID flag MC
extern int    *hgParent; // Hadron origin (not reliable!)


// Gene block =============================================================
extern int gIsub;              // subprocess index
extern double gShat;           // squared total cm energy.
extern double gThat;           // that Mandelstam variable
extern double gUhat;           // that Mandelstam variable
extern double gKsi;            // Ksi variable: (shat + Q^2) / 2*M*nu
extern double gXp;             // x of parton: x/ksi
extern double gZq;             // z of 1st the outgoing quark
extern double gZq2;            // z of 2nd the outgoing quark
extern double gPhiq;           // phi of 1st outgoing quark
extern double gPhiq2;          // phi of 2nd outgoing quark
extern double gAll;            // analysing power
extern double gXi1;            // momentum fraction of the photon parton
extern double gXi2;            // momentum fraction of the nucleon parton
extern double gScale;          // scale in PYTHIA, see mstp(32) and pari
extern double gqPt1;           // p_T of 1 outgoing quark (isub!=99)
extern double gqPt2;           // p_T of other outgoing quark (isub!=99)
extern double gqTheta1;        // theta of 1 outgoing quark (isub!=99)
extern double gqTheta2;        // theta of other outgoing quark (isub!=99)
extern int gKf;                // flavour of outgoing quark
extern int gKf2;               // flavour of 2nd outgoing quark
extern double gDelpovpnucleon; // (deltap/p) in nucleon
extern double gDelpovpgammin;  // idem for photon, min scenario
extern double gDelpovpgammax;  // idem, max scenario

#ifdef GENE_DEBUG
extern double ggE;           // Energy of gamma* (isub!=99)
extern double gqE1;           // Energy of 1 outgoing quark (isub!=99)
extern double gqE2;           // Energy of other outgoing quark (isub!=99)
extern double ggM;           // Mass of gamma* (isub!=99)
extern double gqM1;           // Mass of 1 outgoing quark (isub!=99)
extern double gqM2;           // Mass of other outgoing quark (isub!=99)
extern double gAll_gen;    // analysing power (calculated using xp, zq and phiq as returned by Lepto)
extern double gXp_gen;     // x parton from the proton (generated by Lepto)
extern double gZq_gen;     // z of 1st the outgoing quark (generated by Lepto)
extern double gPhiq_gen;   // phi of 1st outgoing quark (generated by Lepto)
extern double gPhi_lepton; // azimuthal angle between leptons in gamma-P CMS
extern double gKsi_dbg;
extern double gXp_dbg;
extern double gZq_dbg;
extern double gZq2_dbg;
extern double gPhiq_dbg;
extern double gPhiq2_dbg;
#endif //GENE_DEBUG


void InitHighptTree(TTree *t, int nhadrons, int mcON = 0, int geneON = 0);

#endif //__InitTree__
