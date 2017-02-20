#ifndef __SelectHighpts__
#define __SelectHighpts__
#include "PaEvent.h"


using namespace std;

enum {PHASTHIGHPTDEBUG,  // accept 0 to 100 hadrons (all) 
      PHASTHIGHPTPROD2,    // 2 hadrons or more
      PHASTHIGHPTPROD1};  // 1 hadron or more

enum {  
      PHPTNPART,      // 0  not enough particles
      PHPTNOVTX,      // 1  no vertex 
      PHPTNOBPVTX,    // 2  no good primary vertex
      PHPTNOINPART,   // 3  no incoming particle
      PHPTNOMUPRIM,   // 4  no mu' in the PV
      PHPTMUTRACK,    // 5  no mu track 
      PHPTCHI2MU,     // 6  chi2 of Mu
      PHPTCHARGEMUMUP,// 7  Mu Mu' charge >0 
      PHPTMUPTRACK,   // 8  no mu' track
      PHPTCHI2MUP,    // 9  chi2 of Mu
      PHPTXX0MUP,     // 10 XX0 > 30
      PHPTMUPZFIRST,  // 11 zfirst Mu' > 350
      PHPTMUMOM,      // 12 
      PHPTBADBPMC,    // 13 bad MC beam mometum
      PHPTCROSSCELLS, // 14
      PHPTINTARGET,   // 15
      PHPTW,          // 16
      PHPTQ2,         // 17 q2 too small
      PHPTXBJ,        // 18
      PHPTY,          // 19
      PHPTNHADR,      // 20 not enough hadrons
      PHPTMOM1,       // 21 pT1 > 0.7
      PHPTHADZ1,      // 22 z1 > 0.1
      PHPTMOM2,       // 23 pT2 > 0.4
      PHPTHADZ2,      // 24 z2 > 0.1
      PHPTHADZ12,     // 25 z1+z2 < 0.9
      PHPTOK          // 26  event ok
      //
      //PHPTALTTRIG,  // 12 No trigger condition after trigger mask correction
      //PHPTMAGNETH,  // 17 hadron crossed magnet
      //PHPTACC,      // 18 not in defined acceptance
};

struct HptOptions {
float pt0min;
float pt1min;
float pt2min;
float q2min;
float thetamax;
int   mode;
bool  mcON;
bool  geneON;
bool  pythiaON;
bool  muRecoveryON;
bool  inclusiveON;
bool  targetCutON;
float targetR;
float targetY;
int   year;
};

typedef struct HptOptions HptOpt;

int SelectHighpts(PaEvent& e, string optfile);
int SelectHighpts(PaEvent& e, HptOpt opt);
int SelectHighpts(PaEvent& e, float _pt0min, float _pt1min, float _pt2min, float _q2min, float _thetamax, int _mode = PHASTHIGHPTDEBUG, bool _mcON=0, bool _geneON=0, bool _pythiaON=0, bool _muRecoveryON=0, bool _inclusiveON=0);
int SelectHighptsRun(PaEvent& e);

#endif
