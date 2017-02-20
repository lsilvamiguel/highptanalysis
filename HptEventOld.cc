#include "HptEventOld.h"

using namespace std;

ClassImp(HptEventOld)

HptEventOld::HptEventOld()
        : hNmax(2)
{
	Init();
    Clear();
}

HptEventOld::HptEventOld(int hnmax)
        : hNmax(hnmax)
{
	Init();
    Clear();
}

void HptEventOld::Init()
{
	hChi2 = new double[hNmax];
	hQ = new double[hNmax];
	hY_MM1V = new double[hNmax];
	hX_MM1V = new double[hNmax];
	hY_MM1U = new double[hNmax];
	hX_MM1U = new double[hNmax];
	hY_MM1Y = new double[hNmax];
	hX_MM1Y = new double[hNmax];
	hY_DC1 = new double[hNmax];
	hX_DC1 = new double[hNmax];
	hPhi = new double[hNmax];
	hPt = new double[hNmax];
	hPl = new double[hNmax];
	hZ = new double[hNmax];
	hP = new double[hNmax];
	hE = new double[hNmax];
	hPx = new double[hNmax];
	hPy = new double[hNmax];
	hPz = new double[hNmax];
	hTheta = new double[hNmax];
	hXf = new double[hNmax];
	hZlast = new double[hNmax];
	hNhits = new double[hNmax];
	hTime = new double[hNmax];
	hX0 = new double[hNmax];
	hEc = new double[hNmax];
	hHc1 = new double[hNmax];
	hHc2 = new double[hNmax];
	hEc1 = new double[hNmax];
	hEc2 = new double[hNmax];
	hXCClus = new double[hNmax];
	hYCClus = new double[hNmax];
	hZCClus = new double[hNmax];
	hSolx1 = new double[hNmax];
	hSoly1 = new double[hNmax];
	hSolx2 = new double[hNmax];
	hSoly2 = new double[hNmax];
	hSolx3 = new double[hNmax];
	hSoly3 = new double[hNmax];


	// MC block
	hmcP=new double[hNmax];     // momentum    
	hmcTheta=new double[hNmax]; // polar angle w/r gamma   
	hmcPt=new double[hNmax];
	hmcZ=new double[hNmax];
	hmcXf=new double[hNmax];
	hmcPid=new int[hNmax];

	// hgParent=new int[hNmax]; 
}

void HptEventOld::Clear(Option_t* option)
{
    eRun = eSpill = eEvinspill = eTrig = eTrigAlt = vN =  -1000; //int
	eEvNr = 0; //ulong64
    hMinv = vX = vY = vZ = vChi2 = -1000; //double
    kW2 = kQ2 = kXbj = kEgamma = kY = eImu = ePolu = ePold = ePolc = eD = eNN = -1000 ; //double
	bE = bP = bPx = bPy = mE = mP = mPx = mPy = mZlast = mNhits = mTime = mX0 = -1000; //double
	for(int i=0; i < hNmax; i++) {
		hChi2[i] = hQ[i] = hY_MM1V[i] = hX_MM1V[i] = hY_MM1U[i] = hX_MM1U[i] = hY_MM1Y[i] = hX_MM1Y[i] = -1000;
		hY_DC1[i] = hX_DC1[i] = hPhi[i] = hPt[i] = hPl[i] = hZ[i] = hP[i] = hE[i] = hPx[i] = hPy[i] = hPz[i] = hTheta[i] = -1000;
		hXf[i] = hZlast[i] = hNhits[i] = hTime[i] = hX0[i] = hEc[i] = hHc1[i] = hHc2[i] = hEc1[i] = hSolx1[i] = hSoly1[i] = hSolx2[i] = hSoly2[i] = hSolx3[i] = hSoly3[i] =-1000;
	}

	// MC block
    vmcN = -1000; //int
    vmcX = vmcY = vmcZ = -1000; //double
    kmcW2 = kmcQ2 = kmcXbj = kmcEgamma = kmcY = -1000 ; //double
	bmcP = mmcP = -1000; //double
	for(int i=0; i < hNmax; i++) {
		hmcPt[i] = hmcZ[i] = hmcP[i] = hmcTheta[i] = -1000;
		hmcXf[i] = hmcPid[i] = -1000;
		// hgParent[i] = -1000; 
	}
	gIsub = gKf = gKf2 = -1000;
	gShat = gThat = gUhat = gKsi = gXp = gZq = gZq2= gPhiq = gPhiq2 = gAll = gXi1 = gXi2 = gScale = gqPt1 = gqPt2 = -1000;
	gqTheta1 = gqTheta2 = gDelpovpnucleon = gDelpovpgammin = gDelpovpgammax = -1000;
	gAll_gen = gXp_gen = gZq_gen = gPhiq_gen = -1000;

#ifdef GENE_DEBUG
	gPhi_lepton = ggM = gqM1 = gqM2 = ggE = gqE1 = gqE2 = gKsi_dbg = gXp_dbg = gZq_dbg = -1000;
	gZq2_dbg = gPhiq_dbg = gPhiq2_dbg = -1000;
#endif //GENE_DEBUG
}

HptEventOld::~HptEventOld()
{
    Clear();

	if(hChi2) delete [] hChi2;
	if(hQ) delete [] hQ;
	if(hY_MM1V) delete [] hY_MM1V;
	if(hX_MM1V) delete [] hX_MM1V;
	if(hY_MM1U) delete [] hY_MM1U;
	if(hX_MM1U) delete [] hX_MM1U;
	if(hY_MM1Y) delete [] hY_MM1Y;
	if(hX_MM1Y) delete [] hX_MM1Y;
	if(hY_DC1) delete [] hY_DC1;
	if(hX_DC1) delete [] hX_DC1;
	if(hPhi) delete [] hPhi;
	if(hPt) delete [] hPt;
	if(hPl) delete [] hPl;
	if(hZ) delete [] hZ;
	if(hP) delete [] hP;
	if(hE) delete [] hE;
	if(hPx) delete [] hPx;
	if(hPy) delete [] hPy;
	if(hPz) delete [] hPz;
	if(hTheta) delete [] hTheta;
	if(hXf) delete [] hXf;
	if(hZlast) delete [] hZlast;
	if(hNhits) delete [] hNhits;
	if(hTime) delete [] hTime;
	if(hX0) delete [] hX0;
	if(hEc) delete [] hEc;
	if(hHc1) delete [] hHc1;
	if(hHc2) delete [] hHc2;
	if(hEc1) delete [] hEc1;
	if(hEc2) delete [] hEc2;
	if(hXCClus) delete [] hXCClus;
	if(hYCClus) delete [] hYCClus;
	if(hZCClus) delete [] hZCClus;
	if(hSolx1) delete [] hSolx1;
	if(hSoly1) delete [] hSoly1;
	if(hSolx2) delete [] hSolx2;
	if(hSoly2) delete [] hSoly2;
	if(hSolx3) delete [] hSolx3;
	if(hSoly3) delete [] hSoly3;


	// MC block
	if(hmcP) delete [] hmcP;
	if(hmcTheta) delete [] hmcTheta;
	if(hmcPt) delete [] hmcPt;
	if(hmcZ) delete [] hmcZ;
	if(hmcXf) delete [] hmcXf;
	if(hmcPid) delete [] hmcPid;

	// if(hgParent) delete [] hgParent; 

}
