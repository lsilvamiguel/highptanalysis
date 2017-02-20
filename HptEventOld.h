#ifndef HPT_EVENT_OLD_H
#define HPT_EVENT_OLD_H
// #define GENE_DEBUG 

#include "TObject.h"
#include "TLorentzVector.h"

class HptEventOld : public TObject
{
	public:
		HptEventOld();
		HptEventOld(int hnmax);
		~HptEventOld();

		void Clear(Option_t* option="");

		// event information =====================================================
		int       eRun;     // run number
		int       eSpill;   // spill number
		int       eEvinspill;   // event number in spill
		ULong64_t eEvNr;	// Unique event number (for MC - event nr from lepto)
		int       eTrig;    // trigger mask
		int       eTrigAlt;    // corrected trigger mask
		double  eImu;   // incident muon intensity.
		double  ePolu;   // polar upstream
		double  ePolc;   // polar central 
		double  ePold;   // polar downstream 
		double  eD;      // depolarization factor
		double  eNN;    //Neural Network output


		// vertex ================================================================
		double vX;    // X coordinate of vertex
		double vY;    // Y coordinate of vertex
		double vZ;    // Z coordinate of vertex
		int   vN;    // Number of tracks in vertex
		double vChi2; // Chi2 of primary vertex

		// kinematics ============================================================
		double kEgamma;
		double kQ2;     // mass of virtual photon
		double kXbj;    // bjorken variable
		double kY;      // mu relative energy loss  
		double kW2;     // invariant mass of photon nucleon system

		// mu,mup, virtual photon ================================================
		double        bP;
		double        bPx;
		double        bPy;
		double        bE;
		double        mP;
		double        mPx;
		double        mPy;
		double        mE;
		double        mZlast;
		double        mNhits;
		double        mTime;
		double        mX0;

		// hadrons ===============================================================
		int hN;          // number of hadrons
		int hNmax;       // number of hadrons to store
		double hMinv; // Invariant mass of two hadrons
		double*        hChi2; //[hNmax]
		double*        hQ; //[hNmax]
		double*        hY_MM1V; //[hNmax]
		double*        hX_MM1V; //[hNmax]
		double*        hY_MM1U; //[hNmax]
		double*        hX_MM1U; //[hNmax]
		double*        hY_MM1Y; //[hNmax]
		double*        hX_MM1Y; //[hNmax]
		double*        hY_DC1; //[hNmax]
		double*        hX_DC1; //[hNmax]
		double*        hPhi; //[hNmax]
		double*        hPt; //[hNmax]
		double*        hPl; //[hNmax]
		double*        hZ; //[hNmax]
		double*        hP; //[hNmax]
		double*        hE; //[hNmax]
		double*        hPx; //[hNmax]
		double*        hPy; //[hNmax]
		double*        hPz; //[hNmax]
		double*        hTheta; //[hNmax]
		double*        hXf; //[hNmax]
		double*        hZlast; //[hNmax]
		double*        hNhits; //[hNmax]
		double*        hTime; //[hNmax]
		double*        hX0; //[hNmax]
		double*        hEc; //[hNmax]
		double*        hHc1; //[hNmax]
		double*        hHc2; //[hNmax]
		double*        hEc1; //[hNmax]
		double*        hEc2; //[hNmax]
		double*        hXCClus; //[hNmax]
		double*        hYCClus; //[hNmax]
		double*        hZCClus; //[hNmax]
		double*        hSolx1; //[hNmax]
		double*        hSoly1; //[hNmax]
		double*        hSolx2; //[hNmax]
		double*        hSoly2; //[hNmax]
		double*        hSolx3; //[hNmax]
		double*        hSoly3; //[hNmax]


		// MC block
		double vmcX;    // X coordinate of MC vertex
		double vmcY;    // Y coordinate of MC vertex
		double vmcZ;    // Z coordinate of MC vertex
		int    vmcN;      // Number of tracks in MC vertex
		double kmcEgamma; // energy of virtual photon MC
		double kmcQ2;     //  - mass of virtual photon MC
		double kmcXbj;    // bjorken variable MC
		double kmcY;      // mu relative energy loss MC  
		double kmcW2;     // Invariant mass of proton - foton system
		double bmcP;
		double mmcP;
		double *hmcP; //[hNmax]    
		double *hmcTheta; //[hNmax]   
		double *hmcPt; //[hNmax]
		double *hmcZ; //[hNmax]
		double *hmcXf; //[hNmax]
		int    *hmcPid; //[hNmax]

		// Gene block =============================================================
		int gIsub;             // subprocess index
		double gShat;             // squared total cm energy: Q^2 * (1-x_p) / x_p
		double gThat;             // squared total cm energy.
		double gUhat;             // squared total cm energy.
		double gKsi;           // Ksi variable: (shat + Q^2) / 2*M*nu
		double gXp;            // x second parton: x/ksi
		double gZq;            // z first parton
		double gZq2;            // z second parton
		double gPhiq;            // phi of 1st quark
		double gPhiq2;            // phi of second quark
		double gAll;            // analysing power
		double gXi1;            // momentum fraction of the photon parton
		double gXi2;            // momentum fraction of the nucleon parton
		double gScale;          // scale used in PYTHIA
		double gqPt1;
		double gqPt2;
		double gqTheta1;
		double gqTheta2;
		int gKf;            // flavour of outgoing quark
		int gKf2;            // flavour of 2nd outgoing quark
		double gDelpovpnucleon;       // (deltap)/p of the parton coming from nucleon
		double gDelpovpgammin;       // id, for parton coming from gamma (min scenario)
		double gDelpovpgammax;       // id, for parton coming from gamma (max scenario)
		double gAll_gen;    // analysing power (calculated using xp, zq and phiq as returned by Lepto)
		double gXp_gen;     // x parton from the proton (generated by Lepto)
		double gZq_gen;     // z of 1st the outgoing quark (generated by Lepto)
		double gPhiq_gen;   // phi of 1st outgoing quark (generated by Lepto)
		// int    *hgParent; //[hNmax] 


		// Debug

#ifdef GENE_DEBUG
		// Generator block
		double gPhi_lepton; // azimuthal angle between leptons in gamma-P CMS
		double ggM; // mass of virtual photon
		double gqM1; // mass of first outgoing quark
		double gqM2; // mass of second outgoing quark
		double ggE; // energy of virtual photon
		double gqE1; // energy of first outgoing quark
		double gqE2; // energy of second outgoing quark
		double gKsi_dbg;
		double gXp_dbg;
		double gZq_dbg;
		double gZq2_dbg;
		double gPhiq_dbg;
		double gPhiq2_dbg;
#endif //GENE_DEBUG

	private:
		void Init();

		ClassDef(HptEventOld,8) //High pT event
};

#endif

