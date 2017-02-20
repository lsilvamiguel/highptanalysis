#include <iostream>
#include <cmath>
#include <stdio.h>
#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"

void RescaleMom(PaEvent & e, bool faster_mode)
{
// Rescaling of the particle momentum should be done using info directly from 
// NRM measurement (SM2) and eventually Hall probe (SM1)
// Unfortunately this info is not available yet therefore
// in this function average scaling factor was used.
//
//
// It is global function, should be used once per event
// just after begining of UserEventxx.
// Should work both muon and hadron data. 
//
//                 (by Marcin.STOLARSKI@cern.ch)
//
  static bool if_first=1;
  
  static float x_sm2=0;
  static float y_sm2=0;
  static float z_sm2=0;
  static float z_sm1=0;
  static float z_tar=0; 
  static float detdatcorr=1.0000;
  static float sm1corr=1.0000;
  static float sm2corr=1.0000;
  static float sm1sm2corr=1.0000;

  
  if(if_first)
  {
    const PaSetup& setup = PaSetup::Ref();
    PaMagInfo* m = PaSetup::Ref().PtrMagField()->getMagInfo();

    const int MagnetsNumber=PaSetup::Ref().PtrMagField()->getNumOfMags();
  
    x_sm2 = m[MagnetsNumber-1].xcm/10.; // SM2
    y_sm2 = m[MagnetsNumber-1].ycm/10.; // SM2
    z_sm2 = m[MagnetsNumber-1].zcm/10.; // SM2
    z_sm1 = m[MagnetsNumber-2].zcm/10.; // SM1
    z_tar = setup.TargetCenterZ();      // target 


// Field from the field maps @ NMR position (By component)
// 4000A - 1.57904T
// 5000A - 1.91221T
// 2000A - ?.?????T

    float xx=x_sm2-67.5;
    float yy=y_sm2-45.0;
    float zz=z_sm2+1.5;
    float bx,by,bz;
    setup.MagField(xx,yy,zz,bx,by,bz);


    int usedfield=0;
    if(fabs(by)>1.75) usedfield=5000;
    if(fabs(by)>1.45&&fabs(by)<1.75) usedfield=4000;





    if(usedfield==5000&&fabs(fabs(by/1.91221)-1)>0.003) return; // corr done in det.dat?
    if(usedfield==4000&&fabs(fabs(by/1.57905)-1)>0.003) return; //corr done in det.dat?

// - >0.003 in the return statement was used as 2008 hadron data
//  where produced with correction factor 0.998.
//  This was due to the fact that date were collected with current
//  4990A instead of 5000A. Now knowing better the field we should correct
//  these data too.



    switch(usedfield)
    {
      case 4000:
        detdatcorr=fabs(by/1.57905);
        sm2corr=1.0088/detdatcorr;
        sm1sm2corr=1.0088/detdatcorr;
      break;
      case 5000: 
        detdatcorr=fabs(by/1.91221);
        sm2corr=1.0040/detdatcorr;
        sm1sm2corr=1.0040/detdatcorr;
      break;
    }
    if(faster_mode) if_first=0;
  }

// -Scaling factor 1.0088 comes from average of Eric measuremnts
//  @NMR point he obtained 1.0095.
// -detdatcorr should be around 1.000 for data 2002-2006 and
//  0.998 for 1st pre-production of 2008 hadron data.
// -Scaling factor 1.0000 for SM1 is used as we cannot confirm
//  Eric measurement of the field there. 
//  He found out correction factor 1.0045.
//  It looks like our alignment stability in LAS gives additional
//  possible correction factor  up to 1%  (from K0 mass peak) 
// -Tracks which momentum is determined by both SM1 and SM2:
//  for 160 GeV particle the following empirical formula for the
//  correction factor was found:
//  sm1sm2corr= 0.015*sm1corr+ 0.985*sm2corr
//  It was decided to use sm1sm2corr=sm2corr as the impact of sm1corr
//  is very small and sm1corr itself is not a solid number. 



    vector<PaTrack>& TrVec=  e.vTrack(); 
    vector<PaParticle>& PaVec=  e.vParticle();
 
    for(unsigned int i=0; i<TrVec.size();i++)
    {
  
      int TrackType=-1;
      if(TrVec[i].ZLast()<z_tar) TrackType=0;  //beam
      if(TrackType!=0 && TrVec[i].ZFirst()<z_sm1 && TrVec[i].ZLast()<z_sm2) TrackType=1;  //sm1
      if(TrackType!=0 && TrVec[i].ZFirst()>z_sm1 && TrVec[i].ZLast()>z_sm2) TrackType=2;  //sm2
      if(TrackType!=0 && TrVec[i].ZFirst()<z_sm1 && TrVec[i].ZLast()>z_sm2) TrackType=3;  //sm1+sm2


      if(!TrackType) continue; // do nothing with beam track

      vector<PaTPar>& PaTTrVec=TrVec[i].vTPar();
      for(unsigned int j=0; j<PaTTrVec.size();j++)
      { 
        switch(TrackType)
        {
          case 1:
            PaTTrVec[j](5)/=sm1corr;
            break;
          case 2:
            PaTTrVec[j](5)/=sm2corr;
            break;
          case 3:     
            PaTTrVec[j](5)/=sm1sm2corr;
            break;
          // I don't touch covariance matrix. 
        }
      }
    
      if(TrVec[i].iParticle()<0) continue; //some tracks don't have corresponding particle
      vector<PaTPar>& PaTPaVec=PaVec[TrVec[i].iParticle()].vFitPar();   
      for(unsigned int j=0; j<PaTPaVec.size();j++)
      {
        switch(TrackType)
        {
          case 1:
            PaTPaVec[j](5)/=sm1corr;
            break;
          case 2:
            PaTPaVec[j](5)/=sm2corr;
            break;
          case 3:
            PaTPaVec[j](5)/=sm1sm2corr;
            break;
        }
      }
    }

}
