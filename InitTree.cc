#include "InitTree.h"
#include "TTree.h"
#include <iostream>


using namespace std;


// event information =====================================================
int   eRun;
int   eSpill;
int   eEvinspill;
int   eTrig;
double  eImu;
double  ePolu;
double  ePold;
double  eD;


// vertex ================================================================
double vX;
double vY;
double vZ;
int   vN;
double vChi2;


// MC block
double vmcX;
double vmcY;
double vmcZ;
int    vmcN;



// kinematics ============================================================
double kEgamma;
double kQ2;
double kXbj;
double kY;
double kW2;

// MC block
double kmcEgamma;
double kmcQ2;
double kmcXbj;
double kmcY;
double kmcW2;




// beam,mup ===============================================================
double bP;
double bPx;
double bPy;
double bE;
double mP;
double mPx;
double mPy;
double mE;

// MC block
double bmcP;
double mmcP;

double mZlast;
double mNhits;
double mTime;
double mX0;




// hadrons ================================================================
int hN;
double hMinv;
double *hChi2;
double *hQ;
double *hPhi;
double *hPt;
double *hPl;
double *hX_MM1V;
double *hY_MM1V;
double *hX_MM1U;
double *hY_MM1U;
double *hX_MM1Y;
double *hY_MM1Y;
double *hX_DC1;
double *hY_DC1;


double *hZ;

// lab
double *hP;
double *hE;
double *hPx;
double *hPy;
double *hPz;
double *hTheta;

// cms
double *hXf;

// mu id
double *hZlast;
double *hNhits;
double *hTime;
double *hX0;

// calorimeter information
double *hEc;
double *hHc1;
double *hHc2;
double *hEc1;
//L.S.
double *hEc2;  
double *hXCClus;
double *hYCClus;
double *hZCClus;


// MC block
double *hmcP;
double *hmcTheta;
double *hmcPt;
double *hmcZ;
double *hmcXf;
int    *hmcPid;
int    *hgParent;

// Gene block =============================================================
int gIsub;
double gShat;
double gThat;
double gUhat;
double gKsi;
double gXp;
double gZq;
double gZq2;
double gPhiq;
double gPhiq2;
double gAll;
double gXi1;
double gXi2;
double gScale;
double gqPt1;
double gqPt2;
double gqTheta1;
double gqTheta2;
int gKf;
int gKf2;
double gDelpovpnucleon;
double gDelpovpgammin;
double gDelpovpgammax;


// gene block 2 (debug)


#ifdef GENE_DEBUG
double ggE;
double gqE1;
double gqE2;
double ggM;
double gqM1;
double gqM2;
double gAll_gen;
double gXp_gen;
double gZq_gen;
double gPhiq_gen;
double gPhi_lepton;
double gKsi_dbg;
double gXp_dbg;
double gZq_dbg;
double gZq2_dbg;
double gPhiq_dbg;
double gPhiq2_dbg;
#endif

void InitHighptTree(TTree *t, int nhadrons, int mcON, int geneON) {
  
  cout<<"tree initialization. nhadrons : "<<nhadrons<<endl;

  //  t = new TTree("USR1","eff_data_mc");

  hChi2=new double[nhadrons];
  hQ=new double[nhadrons];
  hX_MM1V=new double[nhadrons];
  hY_MM1V=new double[nhadrons];
  hX_MM1U=new double[nhadrons];
  hY_MM1U=new double[nhadrons];
  hX_MM1Y=new double[nhadrons];
  hY_MM1Y=new double[nhadrons];
  hX_DC1=new double[nhadrons];
  hY_DC1=new double[nhadrons];

  hPhi=new double[nhadrons];
  hPt=new double[nhadrons];
  hPl=new double[nhadrons];
  hZ=new double[nhadrons];
  hP=new double[nhadrons];
  hE=new double[nhadrons];
  hPx=new double[nhadrons];
  hPy=new double[nhadrons];
  hPz=new double[nhadrons];
  hTheta=new double[nhadrons];
  hXf=new double[nhadrons];
  hZlast=new double[nhadrons];
  hNhits=new double[nhadrons];
  hTime=new double[nhadrons];
  hX0=new double[nhadrons];
  hEc=new double[nhadrons];
  hHc1=new double[nhadrons];
  hHc2=new double[nhadrons];
  hEc1=new double[nhadrons];
  //L.S.
  hEc2=new double[nhadrons];  
  hXCClus=new double[nhadrons];
  hYCClus=new double[nhadrons];
  hZCClus=new double[nhadrons];


  if(mcON) {
    hmcP=new double[nhadrons];
    hmcTheta=new double[nhadrons];
    hmcPt=new double[nhadrons];
    hmcZ=new double[nhadrons];
    hmcXf=new double[nhadrons];
    hmcPid=new int[nhadrons];
  }
  if(geneON) {
    hgParent=new int[nhadrons];
  }

  t->Branch("eRun", &eRun, "eRun/I");
  t->Branch("eSpill", &eSpill, "eSpill/I");
  t->Branch("eEvinspill", &eEvinspill, "eEvinspill/I");
  t->Branch("eTrig", &eTrig, "eTrig/I");

  t->Branch("eImu", &eImu, "eImu/D");
  t->Branch("ePolu", &ePolu, "ePolu/D");
  t->Branch("ePold", &ePold, "ePold/D");
  t->Branch("eD", &eD, "eD/D");

  t->Branch("vX",  &vX,  "vX/D");
  t->Branch("vY",  &vY,  "vY/D");
  t->Branch("vZ",  &vZ,  "vZ/D");
  t->Branch("vN",  &vN,  "vN/I");
  t->Branch("vChi2",&vChi2,"vChi2/D");

  t->Branch("kEgamma",&kEgamma,"kEgamma/D");
  t->Branch("kQ2",&kQ2,"kQ2/D");
  t->Branch("kXbj",&kXbj,"kXbj/D");
  t->Branch("kY",&kY,"kY/D");
  t->Branch("kW2",&kW2,"kW2/D");

  t->Branch("bP",&bP,"bP/D");
  t->Branch("bPx",&bPx,"bPx/D");
  t->Branch("bPy",&bPy,"bPy/D");
  t->Branch("bE",&bE,"bE/D");
  t->Branch("mP",&mP,"mP/D");  
  t->Branch("mPx",&mPx,"mPx/D");  
  t->Branch("mPy",&mPy,"mPy/D");  
  t->Branch("mE",&mE,"mE/D");

  t->Branch("mZlast",&mZlast,"mZlast/D");
  t->Branch("mNhits",&mNhits,"mNhits/D");
  t->Branch("mTime",&mTime,"mTime/D");
  t->Branch("mX0",&mX0,"mX0/D");

  t->Branch("hN", &hN, "hN/I");
  t->Branch("hMinv", &hMinv, "hMinv/D");
  t->Branch("hChi2",hChi2,"hChi2[hN]/D");
  t->Branch("hQ",hQ,"hQ[hN]/D");
  t->Branch("hX_MM1V",hX_MM1V,"hX_MM1V[hN]/D");
  t->Branch("hY_MM1V",hY_MM1V,"hY_MM1V[hN]/D");
  t->Branch("hX_MM1U",hX_MM1U,"hX_MM1U[hN]/D");
  t->Branch("hY_MM1U",hY_MM1U,"hY_MM1U[hN]/D");
  t->Branch("hX_MM1Y",hX_MM1Y,"hX_MM1Y[hN]/D");
  t->Branch("hY_MM1Y",hY_MM1Y,"hY_MM1Y[hN]/D");
  t->Branch("hX_DC1",hX_DC1,"hX_DC1[hN]/D");
  t->Branch("hY_DC1",hY_DC1,"hY_DC1[hN]/D");

  t->Branch("hPhi",hPhi,"hPhi[hN]/D");
  t->Branch("hPt",hPt,"hPt[hN]/D");
  t->Branch("hPl",hPl,"hPl[hN]/D");
  t->Branch("hZ",hZ,"hZ[hN]/D");
   
  t->Branch("hP",hP,"hP[hN]/D");
  t->Branch("hE",hE,"hE[hN]/D");
  t->Branch("hPx",hPx,"hPx[hN]/D");
  t->Branch("hPy",hPy,"hPy[hN]/D");
  t->Branch("hPz",hPz,"hPz[hN]/D");
  t->Branch("hTheta",hTheta,"hTheta[hN]/D");
    
  t->Branch("hXf",hXf,"hXf[hN]/D");
   
  t->Branch("hZlast",hZlast,"hZlast[hN]/D");
  t->Branch("hNhits",hNhits,"hNhits[hN]/D");
  t->Branch("hTime",hTime,"hTime[hN]/D");
  t->Branch("hX0",hX0,"hX0[hN]/D");
    
  t->Branch("hEc",hEc,"hEc[hN]/D");
  t->Branch("hHc1",hHc1,"hHc1[hN]/D");
  t->Branch("hHc2",hHc2,"hHc2[hN]/D");
  t->Branch("hEc1",hEc1,"hEc1[hN]/D");
  //L.S.
  t->Branch("hEc2",hEc2,"hEc2[hN]/D");
  t->Branch("hXCClus",hXCClus,"hXCClus[hN]/D");
  t->Branch("hYCClus",hYCClus,"hYCClus[hN]/D");
  t->Branch("hZCClus",hZCClus,"hZCClus[hN]/D");

  if(mcON) {
    t->Branch("vmcX",  &vmcX,  "vmcX/D");
    t->Branch("vmcY",  &vmcY,  "vmcY/D");
    t->Branch("vmcZ",  &vmcZ,  "vmcZ/D");
    t->Branch("vmcN",  &vmcN,  "vmcN/I");
    
    t->Branch("kmcEgamma",&kmcEgamma,"kmcEgamma/D");
    t->Branch("kmcQ2",&kmcQ2,"kmcQ2/D");
    t->Branch("kmcXbj",&kmcXbj,"kmcXbj/D");
    t->Branch("kmcY",&kmcY,"kmcY/D");    
    t->Branch("kmcW2",&kmcW2,"kmcW2/D");    

    t->Branch("bmcP",&bmcP,"bmcP/D");
    t->Branch("mmcP",&mmcP,"mmcP/D");  

    t->Branch("hmcP",hmcP,"hmcP[hN]/D");
    t->Branch("hmcTheta",hmcTheta,"hmcTheta[hN]/D");
    t->Branch("hmcPt",hmcPt,"hmcPt[hN]/D");
    t->Branch("hmcZ",hmcZ,"hmcZ[hN]/D");
    t->Branch("hmcXf",hmcXf,"hmcXf[hN]/D");
    t->Branch("hmcPid",hmcPid,"hmcPid[hN]/I");
  }

  if(geneON) {
    t->Branch("gIsub",&gIsub,"gIsub/I");
    t->Branch("gShat",&gShat,"gShat/D");
    t->Branch("gThat",&gThat,"gThat/D");
    t->Branch("gUhat",&gUhat,"gUhat/D");
    t->Branch("gKsi",&gKsi,"gKsi/D");
    t->Branch("gXp",&gXp,"gXp/D");
    t->Branch("gZq",&gZq ,"gZq/D");
    t->Branch("gZq2",&gZq2 ,"gZq2/D");
    t->Branch("gPhiq",&gPhiq ,"gPhiq/D");
    t->Branch("gPhiq2",&gPhiq2 ,"gPhiq2/D");
    t->Branch("gAll",&gAll ,"gAll/D");
    t->Branch("gXi1",&gXi1,"gXi1/D");
    t->Branch("gXi2",&gXi2 ,"gXi2/D");
    t->Branch("gScale",&gScale ,"gScale/D");
    t->Branch("gqPt1",&gqPt1 ,"gqPt1/D");
    t->Branch("gqPt2",&gqPt2 ,"gqPt2/D");
    t->Branch("gqTheta1",&gqTheta1 ,"gqTheta1/D");
    t->Branch("gqTheta2",&gqTheta2 ,"gqTheta2/D");
    t->Branch("gKf",&gKf ,"gKf/I");
    t->Branch("gKf2",&gKf2 ,"gKf2/I");
    t->Branch("gDelpovpnucleon",&gDelpovpnucleon ,"gDelpovpnucleon/D");
    t->Branch("gDelpovpgammin",&gDelpovpgammin ,"gDelpovpgammin/D");
    t->Branch("gDelpovpgammax",&gDelpovpgammax ,"gDelpovpgammax/D");
  
#ifdef GENE_DEBUG
    t->Branch("ggE",&ggE ,"ggE/D");
    t->Branch("gqE1",&gqE1 ,"gqE1/D");
    t->Branch("gqE2",&gqE2 ,"gqE2/D");
    t->Branch("ggM",&ggM ,"ggM/D");
    t->Branch("gqM1",&gqM1 ,"gqM1/D");
    t->Branch("gqM2",&gqM2 ,"gqM2/D");
    t->Branch("gAll_gen",&gAll_gen ,"gAll_gen/D");
    t->Branch("gXp_gen",&gXp_gen,"gXp_gen/D");
    t->Branch("gZq_gen",&gZq_gen ,"gZq_gen/D");
    t->Branch("gPhiq_gen",&gPhiq_gen ,"gPhiq_gen/D");
    t->Branch("gPhi_lepton",&gPhi_lepton ,"gPhi_lepton/D");
    t->Branch("gKsi_dbg",&gKsi_dbg ,"gKsi_dbg/D");
    t->Branch("gXp_dbg",&gXp_dbg ,"gXp_dbg/D");
    t->Branch("gZq_dbg",&gZq_dbg ,"gZq_dbg/D");
    t->Branch("gZq2_dbg",&gZq2_dbg ,"gZq2_dbg/D");
    t->Branch("gPhiq_dbg",&gPhiq_dbg ,"gPhiq_dbg/D");
    t->Branch("gPhiq2_dbg",&gPhiq2_dbg ,"gPhiq2_dbg/D");
#endif
    t->Branch("hgParent",hgParent,"hgParent[hN]/I");
  }

}

