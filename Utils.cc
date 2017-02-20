#include "Utils.h"
#include "TRandom.h"

bool isVetoFired(int eTrig, Double_t Vp, Double_t Vtot, Double_t Vcalo){
  Bool_t isVp=(gRandom->Rndm()<(Vp)); 
  Bool_t isVtot=(gRandom->Rndm()<(Vtot))||(isVp);
  Bool_t isVcalo=(gRandom->Rndm()<(Vcalo))||(isVtot);
  //             IT  MT           LT           OT            CALO            incMT
  //cout<<eTrig<<'\t';
  eTrig=eTrig & (1 + 2*(!isVtot) + 4*(!isVp) + 8*(!isVtot) + 16*(!isVcalo) +256*(!isVtot));
  //cout<<eTrig<<endl;  

  
  if(eTrig==0){
    //cout<<"Event vetoed"<<endl;
    return (true);
  }
  else {
    return (false);
  }
}//end veto

int  IsCrossHodoscopes(const PaTrack& TrackMu1, int TMsk, bool Hadron)
{
	int AltTrigMsk=0;
	int h1,h2;

	if((TMsk&0x1)==0x1)
	{
		h1=TrackMu1.NHitsFoundInDetect("HI04");
		h2=TrackMu1.NHitsFoundInDetect("HI05");
		if(h1>0 && h2>0 && Hadron) AltTrigMsk+=1;
	}

	if((TMsk&0x2)==0x2)
	{
		h1=TrackMu1.NHitsFoundInDetect("HM04");
		h2=TrackMu1.NHitsFoundInDetect("HM05");
		if(h1>0 && h2>0 && Hadron) AltTrigMsk+=2;
	}
	if((TMsk&0x100)==0x100)
	{
		h1=TrackMu1.NHitsFoundInDetect("HM04");
		h2=TrackMu1.NHitsFoundInDetect("HM05");
		if(h1>0 && h2>0) AltTrigMsk+=256;
	}

	if((TMsk&0x4)==0x4)
	{
		h1=TrackMu1.NHitsFoundInDetect("HL04");
		h2=TrackMu1.NHitsFoundInDetect("HL05");
		if(h1>0 && h2>0 && Hadron) AltTrigMsk+=4;
	}

	if((TMsk&0x8)==0x8)
	{
		h1=TrackMu1.NHitsFoundInDetect("HO03");
		h2=TrackMu1.NHitsFoundInDetect("HO04");
		if(h1>0&&h2>0) AltTrigMsk+=8;
	}

	if((TMsk&0x11f)==0x10 && Hadron)
		AltTrigMsk+=16;

	return AltTrigMsk;
}

double DepolFactor(double x, double y, double q2) {
	double m=G3partMass[5];

	double newDKurek = y * (-(2.*m*m*y*y)/(q2*(1.-x*y))-y+2.) /   (((1.-y)*(1.-y)-2.*m*m*y*y/q2+1.) *   sqrt(1.-(4.*m*m*(1.-x)*x*y*y)/(q2*(1.-x*y)*(1.-x*y)))); // calculation from Kyriluk/Kurek

	return newDKurek;
}

void RemoveBadMiddleTrigger(PaEvent& e, const PaTrack& mup_tr) {
  // Extract information about the Hodoscope planes once and store it  in std::map for fast access
  // If run number changes refresh the information
  static bool first = 1;
  static int run_number = 0;
  static double minX = -9999;
  static double Z = 0;
  if(first || run_number != e.RunNum()) {
    run_number = e.RunNum();
    int idet;
    idet = PaSetup::Ref().iDetector("HM05X1_d");
    if(idet>0) {
      const PaDetect& HM05 =  PaSetup::Ref().Detector(idet);
      minX = HM05.Uorig() + 2.5*HM05.Pitch(); // Remove 3 strips
      Z = HM05.Z();
    } else {
      cerr<<"RemoveBadMiddleTrigger ERROR: Detector HM05Y1_d not found!"<<endl;
      exit(1);
    }
    first = 0;
  }

  if((e.TrigMask()&0x102) == 0) return;

  int Npars = mup_tr.NTPar();
  if(Npars == 0) {
    cerr<<"RemoveBadMiddleTrigger ERROR: track does not have associated helixes!"<<endl;
    return;
  }

  PaTPar partr = mup_tr.vTPar(Npars-1); //Track parameters in last measured point
  PaTPar result;
  bool res = partr.Extrapolate(Z, result, false);
  if(result(1) < minX || !res) {
    int new_mask = e.TrigMask() & (~0x102);
    e.vHeader()[1] = new_mask;
  }
}

