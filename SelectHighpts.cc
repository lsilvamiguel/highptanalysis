#define PHASTextgenLEPTO
// #define GENE_DUMP 
// #define GENE_DUMP2 
// #define USE_MC_VERTEX
// One can set also #define GENE_DEBUG in HptEventOld.h

#include "SelectHighpts.h"
#include "HptEventOld.h"
#include "Utils.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"
#include "PaAlgo.h"


#ifdef PHASTinPYTHIA 
#include "PythiaInterface.h"
#endif

#ifdef PHASTinLEPTO 
#include "LeptoInterface.h"
#endif

#ifdef PHASTextgenLEPTO 
#include "LeptoInterface.h"
#endif

#include "G3part.h"

#define NTRAMAX 1000

HptOpt hpt_opt = {
	0.4, //pt0min
	0.4, //pt1min
	0.0, //pt2min
	0.7, //q2min
	-1.0, //thetamax
	PHASTHIGHPTPROD2, //mode
	false, //mcON
	false, //geneON
	false, //pythiaON
	false, //muRecoveryON
	false, //inclusiveON
	true, //targetCutON
	-9999, //targetR
        -9999, //targetY
	2006 //year
};

//2006
// 	1.3, //targetR
// 	1.3, //targetY



int SelectHighpts(PaEvent& e, string optfile) {
	cerr<<"Reading the options from file is not iplemented yet!!"<<endl;
	cerr<<"Using default opts!!"<<endl;
	return SelectHighptsRun(e);
}

int SelectHighpts(PaEvent& e, HptOpt opt) {
	hpt_opt = opt;
	return SelectHighptsRun(e);
}

int SelectHighpts(PaEvent& e, float _pt0min, float _pt1min, float _pt2min, float _q2min, float _thetamax, int _mode, bool _mcON, bool _geneON, bool _pythiaON, bool _muRecoveryON, bool _inclusiveON) {
	hpt_opt.pt0min       = _pt0min;
	hpt_opt.pt1min       = _pt1min;
	hpt_opt.pt2min       = _pt2min;
	hpt_opt.q2min        = _q2min;
	hpt_opt.thetamax     = _thetamax;
	hpt_opt.mode         = _mode;
	hpt_opt.mcON         = _mcON;
	hpt_opt.geneON       = _geneON;
	hpt_opt.pythiaON     = _pythiaON;
	hpt_opt.muRecoveryON = _muRecoveryON;
	hpt_opt.inclusiveON  = _inclusiveON;

	return SelectHighptsRun(e);
}

int SelectHighptsRun(PaEvent& e) {

  //for inclusive case
  if(hpt_opt.inclusiveON)
    if(hpt_opt.mode != PHASTHIGHPTDEBUG){ 
      cerr<<endl;
      cerr<<"hadron mode : PHASTHIGHPTDEBUG"<<endl;
      cerr<<"all hadrons must be selected"<<endl;
      cerr<<"exiting ..."<<endl;
    }
	static bool first(true);
	static TTree *t(NULL); 
	static HptEventOld* event(NULL);
	static int n=-1;
	n++;

	bool debug=false;
	const double M_mu = G3partMass[5];
	const double M_p = G3partMass[14];
	const double M_pi = G3partMass[8];
	double M_charm = 1.5;

	//static const int nhadronsmax = 100;

	int nhadrons = 0;
	switch(hpt_opt.mode) {
		case PHASTHIGHPTDEBUG:
			nhadrons = 100;
			break;
		case PHASTHIGHPTPROD2:
			nhadrons = 2;
			break;
		case PHASTHIGHPTPROD1:
			nhadrons = 1;
			break;
		default:
			break;
	}
	assert(nhadrons);

	if(first){ // histograms and Ntupes booking block

		t = new TTree("USR1","eff_data_mc");
		// InitHighptTree(t, nhadrons, mcON, geneON); 
		event = new HptEventOld(nhadrons);
		t->Branch("Event","HptEventOld",&event);
		cout<<"MESSAGES:"<<endl;
		cout<<"**** pT1: "<<hpt_opt.pt0min      <<" *****"<<endl;
		cout<<"**** pT2: "<<hpt_opt.pt1min      <<" *****"<<endl;
		cout<<"**** Sum pT2: "<<hpt_opt.pt2min      <<" *****"<<endl;
		cout<<"**** Q2: "<<hpt_opt.q2min       <<" *****"<<endl;
		cout<<"**** Hadron mode: "<<hpt_opt.mode        <<" *****"<<endl;
                cout<<"if '0' accept all hdrons"<<endl;
                cout<<"if '1' 2 hi pt hadrons"<<endl;
                cout<<"if '2' 1 hi pt hadron"<<endl;
		cout<<"**** mcON: "<<hpt_opt.mcON        <<" *****"<<endl;
		cout<<"**** GeneOn: "<<hpt_opt.geneON      <<" *****"<<endl;
		cout<<"**** MuRecoveryON: "<<hpt_opt.muRecoveryON<<" *****"<<endl;
		cout<<"**** "<<"+PID"<<" *****"<<endl;
		cout<<"**** "<<hpt_opt.year<<" *****"<<endl;
		cout<<"**** "<<"+NewMuID"<<" *****"<<endl;
		cout<<"**** inclusiveON : "<<hpt_opt.inclusiveON<<" *****"<<endl;
		cout<<"**** targetCutON : "<<hpt_opt.targetCutON<<" *****"<<endl;
		if(hpt_opt.year==2006){  
		  cout<<"**** "<<"+RescaleMom fucntion"<<" *****"<<endl;
		  cout<<"**** "<<"+RemovingBadMT events"<<" *****"<<endl;
		}
		first = false;
	} // end of histogram booking

	bool accepted = false;

	// first filter 

	//*CUT* N outgoing part > 3
	if(hpt_opt.mode == PHASTHIGHPTPROD2) {
		if(  e.NParticle() < 4 ) return PHPTNPART;  // at least mu,mup,h1,h2
		if(! e.NVertex() ) return PHPTNOVTX;        // at least primary vertex
	}
	if(hpt_opt.mode == PHASTHIGHPTPROD1) {
		if(  e.NParticle() < 3 ) return PHPTNPART;  // at least mu,mup,h1
		if(! e.NVertex() ) return PHPTNOVTX;        // at least primary vertex
	}
	if(hpt_opt.mode == PHASTHIGHPTDEBUG) {
		if(  e.NParticle() < 2 ) return PHPTNPART;  // at least mu,mup
		if(! e.NVertex() ) return PHPTNOVTX;        // at least primary vertex
	}

	if(hpt_opt.mcON) {
		assert(e.NMCvertex()); 
		assert(e.NMCtrack()); 
	}

	// store event information
	event->eRun   =  e.RunNum();
	event->eSpill =  e.SpillNum();
	event->eEvinspill =  e.EvInSpill();
	event->eTrig  =  (e.TrigMask()&0x007fffff);

	//if((event->eTrig&0x71f)!=0x200) return PHPTQ2;

	//static const int nhadronsmax = 100;

	// if(eSpill==87 && eEvinspill==14569) debug=true;

	// Veto Deadtime

	/*
	   if(!e.IsMC()) {
	   Double_t Vp=0.0683637;
	   Double_t Vtot=0.1902009;
	   Double_t Vcalo=0.0710521;
	   if(isVetoFired(eTrig, Vp, Vtot, Vcalo)) { return PHPTVETO; }
	   }
	   */

	// scalers
	if(!e.IsMC()) {
		const vector<unsigned int>& scalers = e.vScaler();

		static unsigned lastnmuinspill = 0;
		static double lasttimeinspill = 0;

		unsigned nmuinspill = 0;
		double timeinspill = 0;
		if(scalers.size() >= 6) { 
			nmuinspill = scalers[0] + scalers[1] + scalers[2] + scalers[3] + scalers[4]+scalers[5];
			timeinspill = e.TimeInSpill();
		}    

		if(nmuinspill < lastnmuinspill) {
			// bos - scalers were reset
			lastnmuinspill = 0;
			lasttimeinspill = 0;
		}

		event->eImu = (nmuinspill - lastnmuinspill) / (timeinspill - lasttimeinspill);

		lastnmuinspill = nmuinspill;
		lasttimeinspill = timeinspill;
	}


	// target polarization 
	event->ePolu = event->ePold = 0;
	const vector<Float_t>& polars = PaSetup::Ref().vTargPolar();
	//  assert(polars.size() == 2);
	if(polars.size() == 2) {
		event->ePolu = polars[0]; 
		event->ePold = polars[1]; 
	}
	if(polars.size() == 3) {
		event->ePolu = polars[0]; 
		event->ePolc = polars[1]; 
		event->ePold = polars[2]; 
	}

	// store info. about primary vertex (if found)
	TLorentzVector lvp(0,0,0,M_p); // proton
	TLorentzVector lvmu;
	TLorentzVector lvmup;  
	TLorentzVector lvq;            // photon
	double          pq;            // lvp.lvq
	TVector3       v3q;            // photon
	TVector3       v3perp2qdiff;       // perp to photon, within diffusion plane

	// store info about primary MC vertex (if found)
	TLorentzVector lvmcmu;
	TLorentzVector lvmcmup;  
	TLorentzVector lvmcq;            // photon
	double          mcpq;            // lvp.lvq
	TVector3       v3mcq;            // photon
	TVector3       v3mcperp2qdiff;   // perp to photon, within diffusion plane

	int iprim = -1;
	int viMuPrim = -1;
	
	
	//*CUT* Select Best PV
	iprim = e.iBestPrimaryVertex();
	if(iprim<0) return PHPTNOBPVTX;
	const PaVertex& v = e.vVertex(iprim);

	//assert(v.IsPrimary());

	//NewMuID
	viMuPrim = v.iMuPrim();
	//viMuPrim = v.iMuPrim(true, true, true, false);
	//viMuPrim = v.iMuPrim(true, true, false);
	
	//*CUT* incoming particle
	int inPart= v.InParticle();
	if(inPart<0) return PHPTNOINPART;

	// if no mu prime recover mu' for CaloTrig
	// CAUTION: I think Seb is no d doing this recovery procedure
	
	int mup_rcv = -1;
	//mIsRec=0;

	//*CUT* Mu'
	//if(viMuPrim < 0  ) { //&& !hpt_opt.muRecoveryON
	if(viMuPrim ==-1) { //&& !hpt_opt.muRecoveryON
	  return PHPTNOMUPRIM; // no scattered mu. Skip.
	} 
	// else if ( v.iMuPrim() == -1 ) {
	// 	mup_rcv = PaAlgo::GetMW1ScatMuon(e);
	// 	//cout << "mup_rcv " << mup_rcv << endl;
	// 	if ( mup_rcv > -1) {
	// 		//mIsRec=1;
	// 		const PaTrack& Mu1t = e.vTrack(mup_rcv);    // it's scattered muon
	// 		const PaParticle& Mu1 = e.vParticle(Mu1t.iParticle());    // it's scattered muon
	// 		int tj = Mu1.iTrack();
	// 		if ( mup_rcv!=tj || Mu1.Q() < 0 ) {
	// 			cout << "something wrong!" << endl;
	// 			cout << "particle id " << mup_rcv << ", track id " << tj << endl;
	// 			cout << "charge " << Mu1.Q() << endl;
	// 		}
	// 		viMuPrim = mup_rcv;
	// 		if (Mu1.Q()<0) return PHPTNOMUPRIM; // wrong charge. Skip.
	// 	}
	// 	else return PHPTNOMUPRIM; // no mu prime and no recovered mu prime
	// }


	const PaParticle& mu = e.vParticle(inPart);
	const PaParticle& mup = e.vParticle(viMuPrim);

	// about the scattered muon
	int itmup =  mup.iTrack();
	int itmu =  mu.iTrack();

	const PaTrack& tmup = e.vTrack()[itmup];
	const PaTrack& tmu = e.vTrack()[itmu];


	if(hpt_opt.year==2006){
	  RemoveBadMiddleTrigger(e, tmup);
	  // 	  if((e.TrigMask()&0x007fffff) == 0) return PHPTNOMUPRIM;
	  //LS if((e.TrigMask()&0x71f)==0)return PHPTNOMUPRIM;
 	}

	//*CUT* Mu track
	if(itmu == -1) {
		//cout<<"no track corresponding to scattered mu"<<endl;
	  return PHPTMUTRACK; 
	} 

	//*CUT* Chi2 Mu < 10
	if (tmu.Chi2tot()/(tmu.NHits()-5.)  > 10.) return       PHPTCHI2MU;

	//*CUT* charge of Mu and Mu' > 0
	if(mu.Q() < 0 || mup.Q() < 0) return PHPTCHARGEMUMUP;

	//*CUT* Mu' track
	if(itmup == -1) {
		//cout<<"no track corresponding to scattered mu"<<endl;
	  return PHPTMUPTRACK; 
	} 

	//*CUT* Chi2 Mu' < 10	
	if (tmup.Chi2tot()/(tmup.NHits()-5.) > 10.) return       PHPTCHI2MUP;

	//*CUT* XX0 > 30
	if (tmup.XX0() < 30.) return       PHPTXX0MUP;

	//*CUT* zFirst < 350
	if (tmup.ZFirst() > 350 ) return       PHPTMUPZFIRST;


	// vertex variables 
	event->vX    = v.Pos(0);
	event->vY    = v.Pos(1);
	event->vZ    = v.Pos(2);
	event->vN    = v.NOutParticles(); // number of tracks in vertex 
	event->vChi2 = v.Chi2();

	//assert(mu.IsBeam());
	//assert(mup.IsMuPrim() || mup_rcv);

	if(mu.NFitPar() == 0 || mup.NFitPar() == 0) {
		//cout<<"mu or mup doesn't have fitted pars"<<endl;
		//LS return PHPTMUQUALITY; 
	}





	const PaTPar& pbeam = mu.ParInVtx(iprim);

	lvmu = pbeam.LzVec(M_mu);
	event->bP = lvmu.Vect().Mag();

	//*CUT*  140 < Momtum Mu < 180
	if(event->bP <= 140 || event->bP >= 180) return PHPTMUMOM;

// 	if(tmu.qP() <= 140 || tmu.qP() >= 180){
// 	  std::cout <<"Spill : "<<event->eSpill<<" Event : " <<event->eEvinspill <<" P : "<<tmu.qP()<<std::endl;
//  	  return PHPTMUMOM;}

	// Target Cuts
	
	if(hpt_opt.targetCutON) {
	  int run;
	  if(!e.IsMC()) run = e.RunNum();
	  else {
	    if(hpt_opt.year==2006) run = -3;
	    if(hpt_opt.year==2007) run = -4;
	    if(hpt_opt.year<2005)  run = -2;
	  }
	  
	  //*CUT* beam crosses the target	
	  if( !PaAlgo::CrossCells(mu.ParInVtx(iprim), run, hpt_opt.targetR, hpt_opt.targetY) ) return PHPTCROSSCELLS;
	  
	  // vertex in the target
	  bool up = false;
	  bool ct = false;
	  bool dw = false;
	  
	  up = PaAlgo::InTarget(mu.ParInVtx(iprim), 'U', run, hpt_opt.targetR, hpt_opt.targetY);
	  dw = PaAlgo::InTarget(mu.ParInVtx(iprim), 'D', run, hpt_opt.targetR, hpt_opt.targetY);
	  if(hpt_opt.year > 2004) ct = PaAlgo::InTarget(mu.ParInVtx(iprim), 'C', run, hpt_opt.targetR, hpt_opt.targetY);
	  
	  //*CUT* Vertex in target
	  if( !( up || dw || ct )) return PHPTINTARGET;
	}

	if( PaSetup::Ref().NDetectors() )
		event->mZlast = tmup.ZLast();
	else 
		event->mZlast = 0;

	event->mNhits = (double) tmup.NHits();
	event->mTime = tmup.MeanTime();
	event->mX0 = tmup.XX0();

	PaTPar partr_first = tmup.vTPar(0);      //Track parameters in first measured point
	PaTPar  pmup_out;

	// calculating inclusive kinematic variables

	event->bPx = lvmu.Vect().Px();
	event->bPy = lvmu.Vect().Py();
	event->bE = lvmu.E();
	lvmup = mup.ParInVtx(iprim).LzVec(M_mu); 
	event->mP = lvmup.Vect().Mag();    
	event->mPx = lvmup.Vect().Px();
	event->mPy = lvmup.Vect().Py();
	event->mE = lvmup.E();

	//if(e.IsMC() && (event->bP > 180 || event->bP < 140)) return PHPTBADBPMC;

	lvq = lvmu - lvmup;
	v3q = lvq.Vect();

	event->kQ2 = -lvq.M2();
	TLorentzVector P_q;
	P_q = lvp + lvq;
	event->kW2 = P_q.M2();
	
	//*CUT* W>5
	if(sqrt(event->kW2) < 5) return PHPTW;
	
	pq = lvp.Dot(lvq);
	double pk = lvp.Dot(lvmu);

	event->kEgamma = lvq.E();
	event->kXbj = event->kQ2/(2*pq);
	event->kY = pq/pk; 
	event->eD = DepolFactor(event->kXbj, event->kY, event->kQ2);
	// v3perp2q = q x (mu x mup)
	v3perp2qdiff = lvmu.Vect().Cross( lvmup.Vect() );
	v3perp2qdiff = v3perp2qdiff.Unit();


// 	if(e.UniqueEvNum() ==  229132378057379LL) {
// 	  cout <<"\t"<< e.UniqueEvNum() <<"\t"<< event->bE - event->mE<<endl<<"\t"<<event->kY<<endl;
// 	}


	if(hpt_opt.mcON) {
#ifdef USE_MC_VERTEX
		const PaMCvertex& vmc = e.vMCvtx(0);

		// first vertex in the list is primary by CG convention
		assert( vmc.IsPrimary() );

		vmcX   = vmc.Pos(0);
		vmcY   = vmc.Pos(1);
		vmcZ   = vmc.Pos(2);
		vmcN   = vmc.NMCtrack() + 1; // +1 for the incoming mu.

		int imcmu = -1;
		int imcmup = -1;

		// to get mu and mup, no other way than looping on mc tracks. 
		// in the MC primary vertex. 
		for(int imct=0; imct<vmc.NMCtrack(); imct++) {
			const PaMCtrack& mct = e.vMCtrk( vmc.iMCtrack(imct) );

			if( mct.Pid() != 5)  continue; // mu or mup

			if( mct.IsPileup() ) continue;
			else if( mct.IsBeam() ) imcmu = imct;
			else imcmup = imct;
		}
#else
		int imcmu = tmu.iMCtrack();
		int imcmup = tmup.iMCtrack();
#endif //USE_MC_VERTEX

		if(imcmu == -1 || imcmup == -1) {
			cout<<"SelectHighpts : warning : missing MC mu or mu prime."<<endl;

			event->kmcEgamma = -1; 
			event->kmcQ2 = -1;     
			event->kmcXbj = -1;    
			event->kmcY = -1;
			event->kmcW2 = -1;
		}
		else {
			const PaMCtrack& mcmu  = e.vMCtrk( imcmu );
			const PaMCtrack& mcmup = e.vMCtrk( imcmup );

			lvmcmu = mcmu.ParInVtx().LzVec(M_mu);
			lvmcmup = mcmup.ParInVtx().LzVec(M_mu); 
			event->bmcP = lvmcmu.Vect().Mag();
			event->mmcP = lvmcmup.Vect().Mag();    
			lvmcq = lvmcmu - lvmcmup;
			v3mcq = lvmcq.Vect();

			event->kmcQ2 = -lvmcq.M2();

			mcpq = lvp.Dot(lvmcq);
			double mcpk = lvp.Dot(lvmcmu);

			event->kmcEgamma = lvmcq.E();
			event->kmcXbj = event->kmcQ2/(2*mcpq);
			event->kmcY = mcpq/mcpk; 
			TLorentzVector mcP_q;
			mcP_q = lvp + lvmcq;
			event->kmcW2 = mcP_q.M2();

			v3mcperp2qdiff = lvmcmu.Vect().Cross( lvmcmup.Vect() );
			v3mcperp2qdiff = v3mcperp2qdiff.Unit();

		}
	}


	// if(event->kXbj>=0 && event->kXbj<=1 &&
	// 		event->kY>=0.1 && event->kY<=0.9 &&
	// 		event->kQ2>=0 ) {

	// 	accepted = true; 
	// 	//cout<< e.UniqueEvNum()<<endl;
	// }

	// -----------------------------------------------------------------

	// something is wrong with the primary vertex

	//	if( !accepted ) return PHPTKINECLEAN;  //
	
	//*CUT* Q^2 > 1
	if(event->kQ2 < hpt_opt.q2min) return PHPTQ2; // q2 too small
	
	//*CUT* 0.003 < x < 0.7
	if(event->kXbj < 0.003 || event->kXbj > 0.7)   return PHPTXBJ;

	//*CUT* 0.1 < y < 0.9
	if(event->kY < 0.1 || event->kY > 0.9) return PHPTY;

	//assert(iprim != -1);

	// -----------------------------------------------------------------
	// sort hadrons in primary vertex by decreasing Pt

	map<double, int, greater<double> > hptids;
	const PaVertex& primv = e.vVertex(iprim);
	if(!hpt_opt.inclusiveON) assert(primv.IsPrimary() );
	if(debug) {
		cout<<"sorting hadrons"<<endl;
		cout<<primv.NOutParticles()<<" particles in primary vertex"<<endl;
	}
	//cout<<"sorting hadrons: "<<primv.NOutParticles()<<endl;
	
	for( int ip=0; ip!=primv.NOutParticles(); ip++) {

		if(debug) cout<<ip<<"\t";

		int pid = primv.iOutParticle(ip);

		if(!hpt_opt.inclusiveON && pid == viMuPrim) continue; // this is the mu'    

		const PaParticle &ptc = e.vParticle() [ pid ];
		if(!hpt_opt.inclusiveON) assert(!ptc.IsBeam());

		//**if(!hpt_opt.inclusiveON && (ptc.PID()==5 || ptc.PID()==6)) continue;

		// from now on, we know that the particle is not the mu'.
		// if another mu' is present in the vertex, event is rejected.
		//LS if(!hpt_opt.inclusiveON &&  ptc.IsMuPrim() ) return PHPT2MUP;

		// nb other particles tagged as mu' are accepted as hadrons.
		// Mon Nov 24 14:30:52 MET 2003 removed the condition hadron is not mu'

		int itrack =  ptc.iTrack();
		const PaTrack& tck = e.vTrack()[itrack];

		if(!hpt_opt.inclusiveON && itrack == -1) {
			if(debug) cout<<"has no track associated"<<endl;
			continue; 
		} 

		if(hpt_opt.mcON) {
			// if mc block is on, do not consider tracks that are not 
			// associated to an MC track
			if(!hpt_opt.inclusiveON && tck.iMCtrack() == -1) 
				continue;
		}

		if(!hpt_opt.inclusiveON && tck.XX0() >= 10) continue;
 		//**if(!hpt_opt.inclusiveON && tck.Chi2tot()/(tck.NHits()-5) > 10) continue;
 		if(!hpt_opt.inclusiveON && tck.Chi2tot()/tck.Ndf() >= 10) continue;
		if(!hpt_opt.inclusiveON && tck.ZFirst() >= 350) continue;
		if(!hpt_opt.inclusiveON && (tck.ZLast() >= 3300 || tck.ZLast() <= 350)) continue;


		const PaTPar& par = ptc.ParInVtx(iprim);      
		//**if(!hpt_opt.inclusiveON) assert(par.HasMom());


		TVector3 v3h = par.Mom3();
		double pt = v3h.Pt(v3q);

		if(debug) cout<<pt<<endl;
		hptids[ pt ] = pid;    
		//cout<<"ID "<<pid<<" Pt: "<<pt<<endl;

	}

	// Prepare corrected trigger mask - check that both trrigger planes fired
	// For seminiclusive case check if we have hadron (in this case secong outgoing part in PV)
	// event->eTrigAlt = IsCrossHodoscopes(tmup, event->eTrig, hptids.size());

	// Check which triggers fired
	//if( event->eTrigAlt==0 ) { return PHPTALTTRIG; }

	// -----------------------------------------------------------------
	// something is wrong with the hadrons
	if(!hpt_opt.inclusiveON && hpt_opt.mode == PHASTHIGHPTDEBUG) {
		if(hptids.size() < 1 && (hpt_opt.pt0min || hpt_opt.pt1min || hpt_opt.pt2min)) return PHPTNHADR; // not enough hadrons
		if(hptids.size()) {

			typedef map<double, int, greater<double> >::iterator IM;   
			map<double, int, greater<double> >::iterator im = hptids.begin();

			double pt2sum = 0;
			int h1=0 ,h2=0; // 1 hadron or more
			int ih = 0;
			for(IM im=hptids.begin(); im!=hptids.end(); im++) {
				//  im++;

				if(ih==0) {if(im->first < hpt_opt.pt0min) return PHPTMOM1;} // mom of 1st hadron too small
				if(ih==1) {if(im->first < hpt_opt.pt1min) return PHPTMOM2;} // mom of 2nd hadron too small
				//     cout << "numero du hadron:" << ih << endl;
				pt2sum += im->first * im->first;
				//    if(++ih == 1 && im++==hptids.end()) {h1 = 1; break;im--;}
				if(ih==1) {h2 = 1; break;} // we consider only the 1st 2 hadrons
				ih++;
			}   
			if (h2==1) {if(pt2sum < hpt_opt.pt2min) return PHPTMOM2;} // 2 hadrons
			if (h2==0) {                                      // 1 hadron
				if(pt2sum < 0.5*hpt_opt.pt2min) return PHPTMOM2;}
			// if(pt2sum < pt2min) return PHPTMOM2; // sum of squared pt too small
		}
	} 


	if(hpt_opt.mode == PHASTHIGHPTPROD2) {
	  
	  //*CUT* N hadrons >= 2
	  if(hptids.size() < 2){
	    //cout<<"size: "<<hptids.size()<<endl;
	    return PHPTNHADR; // not enough hadrons
	  }


		typedef map<double, int, greater<double> >::iterator IM;   
		map<double, int, greater<double> >::iterator im = hptids.begin();

		double pt2sum = 0;
		int ih = 0;

		double z1=0;
		double z2=0;
		  
		for(IM im=hptids.begin(); im!=hptids.end(); im++) {
		           
		  //LS
		  if(ih==0) {
		    //*CUT* pT1 > 0.7
		    //cout<<"** 1 pt: "<<im->first <<endl;
		    if(im->first < hpt_opt.pt0min) 
		      return PHPTMOM1; // mom of 1st hadron too small
		    
		    		    
		    const vector<PaParticle>& parts = e.vParticle();
		    const PaParticle& ptc = parts[im->second];
		    int itrack =  ptc.iTrack();
		    assert(itrack != -1); 
		    
		    const PaTrack& tck = e.vTrack()[itrack];

		    // hadron is assumed to be a pion...
		    TLorentzVector lvh = ptc.ParInVtx(iprim).LzVec(M_pi);
		    
		    z1 = lvp.Dot(lvh) / pq;

		    //cout<<"** 1 z: "<<z1 <<endl;
		    
		    //*CUT* z1 > 0.1
		    if( z1 <= 0.1 ) return PHPTHADZ1;

		  }

		  if(ih==1) {
		    //*CUT* pT1 > 0.4

		    //cout<<"** 2 pt: "<<im->first <<endl;
		    if(im->first < hpt_opt.pt1min) 
		      return PHPTMOM2; // mom of 2nd hadron too small

		    const vector<PaParticle>& parts = e.vParticle();
		    const PaParticle& ptc = parts[im->second];
		    int itrack =  ptc.iTrack();
		    assert(itrack != -1); 
		    
		    const PaTrack& tck = e.vTrack()[itrack];

		    // hadron is assumed to be a pion...
		    TLorentzVector lvh = ptc.ParInVtx(iprim).LzVec(M_pi);

		    z2 = lvp.Dot(lvh) / pq;
		    
		    //cout<<"** 2 z: "<<z1 <<endl;
		    //*CUT* z2 > 0.1
		    if( z2  <= 0.1 ) return PHPTHADZ2;

		  }

		  pt2sum += im->first * im->first;

		  if(++ih == 2) break; // we consider only the 1st 2 hadrons
		}   
		
		//*CUT* z1+z2 < 0.9
		if( z1+z2 >= 0.9 ) return PHPTHADZ12;
		
		//LS if(pt2sum < hpt_opt.pt2min) return PHPTMOM2; // sum of squared pt too small
	} 



	if(hpt_opt.mode == PHASTHIGHPTPROD1) {

		if(hptids.size() < 1) return PHPTNHADR; // not enough hadrons

		typedef map<double, int, greater<double> >::iterator IM;   
		map<double, int, greater<double> >::iterator im = hptids.begin();


		double pt2sum = 0;
		int ih = 0;
		for(IM im=hptids.begin(); im!=hptids.end(); im++) {
			//  im++;

			if(im->first < hpt_opt.pt0min) return PHPTMOM1; // mom of hadrons too small
			pt2sum += im->first * im->first;
			if(++ih == 1) break; // we consider only the 1st hadron !!
		}   

		if(pt2sum < hpt_opt.pt2min) return PHPTMOM2; // sum of squared pt too small

		//    if(hptids.size() < 1) return PHPTNHADR; // not enough hadrons

		//   map<double, int, greater<double> >::iterator im = hptids.begin();
		//  im++;
		//  if(im->first < ptmin) return PHPTMOM; // mom of hadrons too small
	} 

// 	//L.S. to check thresholds
// 	event->eNClus= e.NCaloClus();
		  
// 	int NCHc1=0;
// 	int NCHc2=0;

// 	for(int ical = 0; ical < e.NCaloClus(); ical ++) {
// 	  //cout <<++iclus<<"\t"<<e.NCaloClus()<<endl;
// 	  const PaCaloClus& cclus =  e.vCaloClus(ical);
	  
// 	  if(hpt_opt.year > 2004){
// 	    //For HCAL1
// 	    if(cclus.Z()>1200 && cclus.Z()<1350) NCHc1++;
// 	    //For HCAL2 
// 	    if(cclus.Z()>3450 && cclus.Z()<3800) NCHc2++;
// 	  }else{
// 	    //For HCAL1
// 	    if(cclus.Z()>1150 && cclus.Z()<1400) NCHc1++;
// 	    //For HCAL2 
// 	    if(cclus.Z()>3450 && cclus.Z()<3800) NCHc2++;
// 	  }
// 	  //cout <<Eclus1<< "\t"<<Eclus2<<endl;
// 	}

// 	event->eNClusHc1=NCHc1;
// 	event->eNClusHc2=NCHc2; 

// 	//L.S. to check thresholds
// 	double thres=8.0;
// 	double Eclus1=0;
// 	double Eclus2=0;
// 	//int iclus=0;
// 	for(int ical = 0; ical < e.NCaloClus(); ical ++) {
// 	  //cout <<++iclus<<"\t"<<e.NCaloClus()<<endl;
// 	  const PaCaloClus& cclus =  e.vCaloClus(ical);
// 	  if(hpt_opt.year > 2004){
// 	    //HCAL1
// 	    if(cclus.Z()>1200 && cclus.Z()<1350 && Eclus1<cclus.E())
// 	      Eclus1=cclus.E();
// 	    //HCAL2
// 	    if(cclus.Z()>3450 && cclus.Z()<3800 && Eclus2<cclus.E())
// 	      Eclus2=cclus.E();
// 	  }else{
// 	    //HCAL1
// 	    if(cclus.Z()>1150 && cclus.Z()<1400 && Eclus1<cclus.E())
// 	      Eclus1=cclus.E();
// 	    //HCAL2
// 	    if(cclus.Z()>3450 && cclus.Z()<3800 && Eclus2<cclus.E())
// 	      Eclus2=cclus.E();
// 	  }
// 	  //cout <<Eclus1<< "\t"<<Eclus2<<endl;
// 	}
	
// //


// 	//if(Eclus1<thres && Eclus2<4.5) return PHPTBADH;
// 	if(Eclus1<7.0 && Eclus2<thres) return PHPTBADH;


	// -----------------------------------------------------------------
	// filling hadron block

	int hN = 0;
	event->hMinv = 0;
	TLorentzVector vect;
	double SumZ = -1;

	const vector<PaParticle>& parts = e.vParticle();
	for(map<double,int,greater<double> >::iterator im=hptids.begin();
			im!=hptids.end(); im++) {

		const PaParticle& ptc = parts[im->second];

		// if no track associated to ptcle, continue
		int itrack =  ptc.iTrack();
		assert(itrack != -1); 

		const PaTrack& tck = e.vTrack()[itrack];
		//#define MMstudies
#ifdef MMstudies
		const PaTPar& par = ptc.ParInVtx(iprim);
		//   PaTPar& par1 = par;
		PaTPar  par1;
		par1(0) = 0; par1(1) = 0; par1(2) = 0; par1(3) = 0; par1(4) = 0; par1(5) = 0;

		par.Extrapolate(140.9,par1,true);
		event->hX_MM1V[hN]=par1(1);
		event->hY_MM1V[hN]=par1(2);
		par.Extrapolate(142.3,par1,true);
		event->hX_MM1U[hN]=par1(1);
		event->hY_MM1U[hN]=par1(2);
		par.Extrapolate(151.55,par1,true);
		event->hX_MM1Y[hN]=par1(1);
		event->hY_MM1Y[hN]=par1(2);
		par.Extrapolate(265.,par1,true);
		event->hX_DC1[hN]=par1(1);
		event->hY_DC1[hN]=par1(2);

#endif

		// nb here we know the track has a momentum !
		event->hChi2[hN] = tck.Chi2tot()/(tck.NHits()-5.);

		event->hQ[hN] = ptc.Q();

		if( PaSetup::Ref().NDetectors() )
			event->hZlast[hN] = tck.ZLast();
		else
			event->hZlast[hN] = 0;

		event->hNhits[hN] = (double) tck.NHits();
		event->hTime[hN] = tck.MeanTime();
		event->hX0[hN] = tck.XX0();

		// we also know it's connected to prim vtx
		//cout<<"iprim "<<iprim<<endl;
		TVector3 v3h = ptc.ParInVtx(iprim).Mom3();

		// hadron is assumed to be a pion...
		TLorentzVector lvh = ptc.ParInVtx(iprim).LzVec(M_pi);

		TVector3 v3perp2qhadr = v3q.Cross(v3h);
		v3perp2qhadr = v3perp2qhadr.Unit();

		if( v3q.Dot( v3perp2qdiff.Cross(v3perp2qhadr) ) >= 0)
			event->hPhi[hN] = acos(v3perp2qhadr * v3perp2qdiff);
		else 
			event->hPhi[hN] = -acos(v3perp2qhadr * v3perp2qdiff);

		// in the lab
		event->hP[hN] = v3h.Mag();
		event->hE[hN] = lvh.E();
		event->hPx[hN] = lvh.Px();
		event->hPy[hN] = lvh.Py();
		event->hPz[hN] = lvh.Pz();

		event->hTheta[hN] = v3h.Theta();
		event->hXf[hN] = PaAlgo::Xf(lvmu, lvmup, lvh);
		
		// *** NEW ***


		//para apagar

// 		const PaTPar& par = ptc.ParInVtx(iprim);
// 		PaTPar wind;
// 		wind(0)= wind(1)= wind(2)= wind(3)= wind(4)= wind(5)=0;

// 		//par.Extrapolate(77,wind,true);
// 		par.Extrapolate(118.4,wind,true);

// 		//if(sqrt(wind(1)*wind(1)+wind(2)*wind(2)) > 11 ){ 
// 		if(sqrt(wind(1)*wind(1)+wind(2)*wind(2)) > 14.9 ){ 
// 		  //cout << "out \t";
// 		  continue;
// 		}

		// Pt wrt virtual photon
		event->hPt[hN] = im->first;
		//if(hN<2)cout <<"pT["<<hN<<"]: "<<event->hPt[hN]<<endl;
		event->hPl[hN] = sqrt(event->hP[hN]*event->hP[hN]-event->hPt[hN]*event->hPt[hN]);
		event->hZ[hN] = lvp.Dot(lvh) / pq;

		// calorimeter information (energy deposited in all calorimeters)
		event->hEc[hN] = 0;
		event->hHc1[hN] = 0;
		event->hHc2[hN] = 0;
		event->hEc1[hN] = 0;
		event->hEc2[hN] = 0;

// 		//L.S. to check thresholds
// 		//double thres=5.0;
// 		double Eclus1=-10;
// 		double Eclus2=-10;
// 		//int iclus=0;

// 		event->hNCluster[hN]= ptc.NCalorim();
		for(int ical = 0; ical < ptc.NCalorim(); ical ++) {
		  const PaCaloClus& cclus =  e.vCaloClus()[ ptc.iCalorim(ical) ];
		  //event->hClusSize[hN]= cclus.Size();		  	
		  //L.S.
		  const string& nam = cclus.CalorimName();
 		  if(hpt_opt.year > 2004){
// 		    if(cclus.Z()>1050 && cclus.Z()<1150)   event->hEc1[hN] += cclus.E();
// 		    if(cclus.Z()>1200 && cclus.Z()<1350)   event->hHc1[hN] += cclus.E();
// 		    if(cclus.Z()>3200 && cclus.Z()<3400)   event->hEc2[hN] += cclus.E();
// 		    if(cclus.Z()>3450 && cclus.Z()<3800)   event->hHc2[hN] += cclus.E();
 		    if(nam.find("EC01P1")!=0)   event->hEc1[hN] += cclus.E();
 		    if(nam.find("HC01P1")!=0)   event->hHc1[hN] += cclus.E();
 		    if(nam.find("EC02P1")!=0)   event->hEc2[hN] += cclus.E();
 		    if(nam.find("HC02P1")!=0)   event->hHc2[hN] += cclus.E();
		  }else{
		    if(cclus.Z()>3300 && cclus.Z()<3450)   event->hEc2[hN] += cclus.E();
		    if(cclus.Z()>1150 && cclus.Z()<1400)   event->hHc1[hN] += cclus.E();
		    if(cclus.Z()>3450 && cclus.Z()<3800)   event->hHc2[hN] += cclus.E();
		  }
		  event->hXCClus[hN] = cclus.X();
		  event->hYCClus[hN] = cclus.Y();
		  event->hZCClus[hN] = cclus.Z();
		  
		}
		event->hEc[hN] = event->hHc1[hN]+event->hHc2[hN];
		  //+event->hEc1[hN];
		
		if(! ptc.NCalorim() ){
		  event->hEc[hN] = -1;
		}

// 		if(hN<2 &&
// 		   !(((event->hEc[hN])>0&&(event->hEc[hN])/event->hP[hN]>0.3)||
// 		     (event->hEc[hN]<=0&&event->hZlast[hN]<4000)))
// 		  return PHPTBADH; 		           
		 		//cout<<"IN"<<endl;

		if(hpt_opt.mcON) {
			int imc = tck.iMCtrack();
			assert( imc != -1 );
			const PaMCtrack& mctck = e.vMCtrk(imc);
			int igen = mctck.iGenParticle() - 1;
			// cout<<igen<<'\t'; 
			TLorentzVector lvmch = mctck.ParInVtx().LzVec(M_pi);

			TVector3 v3mcp = lvmch.Vect();


			event->hmcP[hN] = v3mcp.Mag();
			event->hmcTheta[hN] = v3mcp.Theta(); 
			event->hmcPt[hN] = v3mcp.Pt(v3mcq);
			event->hmcZ[hN] = lvp.Dot(lvmch) / mcpq;
			event->hmcXf[hN] = PaAlgo::Xf(lvmcmu, lvmcmup, lvmch);
			event->hmcPid[hN] = mctck.Pid();
			// cout<<mctck.Pid()<<endl;
		}


//  		if(hN<2 && (e.UniqueEvNum() ==  229132378057379LL)) {
//  		  cout<<e.UniqueEvNum()<<"\t"<<hN<<"\t"<<event->hXf[hN]<<"\t"<<event->hZ[hN]<<"\t"<<event->hChi2[hN]<<"\t"<<event->hZlast[hN]<<"\t"<<event->hP[hN]<<endl;
		  
//  		}


		if(!hpt_opt.inclusiveON && hN<2 && (
					event->hXf[hN]<0 || event->hXf[hN]>1 ||
					event->hZ[hN]<0 || event->hZ[hN]>1 || 
					event->hChi2[hN]>20 ||
					(event->hZlast[hN] < 400 ||
					 event->hZlast[hN] > 4000)) ) {


			// cerr<<"bad hadron !"<<endl;
			// this hadron doesn't look good throwing away whole event

			// commented for comparison with Sonja

			//LS return PHPTBADH;
		}
 		// if( hN<2) {
		  
		//   // Remove solenoid cut for 2006
		//   const PaTPar& par1 = tck.vTPar(0);
		//   PaTPar  par_out;

		//   if(hpt_opt.year < 2006){
		//     par1.Extrapolate(118.4,par_out,true);
		//     event->hSolx1[hN] = par_out(1);
		//     event->hSoly1[hN] = par_out(2);
		    
		//     if(!hpt_opt.inclusiveON && sqrt(par_out(1)*par_out(1)+par_out(2)*par_out(2)) > 14.9 ) return PHPTMAGNETH;
		//   }
		//   else {
		//     par1.Extrapolate(130.5,par_out,true);
		//     event->hSolx1[hN] = par_out(1);
		//     event->hSoly1[hN] = par_out(2);
		    
		//     if(!hpt_opt.inclusiveON && sqrt(par_out(1)*par_out(1)+par_out(2)*par_out(2)) > 35 ) return PHPTMAGNETH;
		//   }

 		// }

		// if(hN==0) {
		// 	vect = parts[im->second].ParInVtx(iprim).LzVec(M_pi);
		// 	SumZ = event->hZ[hN];
		// } else if(hN==1) {
		// 	vect += parts[im->second].ParInVtx(iprim).LzVec(M_pi);
		// 	SumZ += event->hZ[hN];
		// 	event->hMinv = vect.M();
		// 	if(!hpt_opt.inclusiveON && SumZ >= 0.95) return PHPTBADH; 
		// }

		hN++;
		if(hN == nhadrons) break;
	}
	event->hN = hN;


	// -----------------------------------------------------------------
	// end of event

// 	if(hpt_opt.thetamax > 0) {
// 		// first 2 hadrons must have theta<thetamax
// 		// magnet acceptance cut

// 		for(int ih=0; ih<hN; ih++) {
// 			if(event->hTheta[ih] > hpt_opt.thetamax)   return PHPTACC;
// 		}
// 	}  

// 	// event has been accepted. extracting some information from generator
// 	// from pythya to calculate a_ll
// 	//TODO test with reading Ahmed lepto
// 	//TODO would be good to use the same events in reading and creating lepto by Ahmed
// 	if(hpt_opt.geneON && e.IsMC() ) { // Generator block

// 		// these are the position of the relevant particules in event record
// 		int igamma = -1;
// 		int imu = -1;
// 		int imup = -1;
// 		int iparton = -1;
// 		int ioutq = -1;
// 		int ioutq2 = -1;

// 		// these are what we want to calculate 
// 		int isub = -1;       // subprocess index
// 		int kftarget = -10;
// 		int kfa = -10;       // flavour of incoming parton (resolved case) from gamma
// 		int kfb = -10;
// 		int kf = -10;         // flavour of outgoing quark
// 		int kf2 = -10;        // flavour of 2nd outgoing quark
// 		double shat = -1;    // energy in com 
// 		double uhat = -1;
// 		double that = -1;

// 		double xp = -1;      // x parton 
// 		double zq = -1;      // E*_q/E*_l
// 		double zq2 = -1;      // E*_q/E*_l
// 		double ksi = -1;
// 		double phi = -10;    // azymuthal angle of outgoing quark w/r leptonic plane
// 		double phi2 = -10;   // azymuth. angle of other quark; not the same(k_T)!
// 		double phi_lepton = -10;
// 		double all = -10;    // a_ll
// 		// double allK = -10;
// 		double xi1 = -1;     // xi parton in photon
// 		double xi2 = -1;     // xi parton in nucleon
// 		double scalePYTHIA = -1; // scale used in PYTHIA
// 		double qPt1;
// 		double qPt2;
// 		double qTheta1;
// 		double qTheta2;
// 		double delpovpnucleon = -10;
// 		double delpovpgammin = -10;
// 		double delpovpgammax = -10;

// 		TLorentzVector lvg;        // gamma
// #ifdef PHASTextgenLEPTO
// 		TLorentzVector lvmu_g;        // generated mu
// 		TLorentzVector lvmup_g;        // generated mu'
// #endif
// 		TLorentzVector lvpfromg;   // parton from gamma
// 		TLorentzVector lvt;        // target
// 		TLorentzVector lvcom;      // center of mass
// 		TLorentzVector lvoutq;     // first outgoing particle (quark or gluon)
// 		TLorentzVector lvoutq2;    // second outgoing particle (quark or gluon)
// 		TLorentzVector lvingamma;
// 		TLorentzVector lvinparton;
// 		TVector3 comboost;         // center of mass boost vector


// #ifdef PHASTinPYTHIA // ----------------------------------------  
// #define INEXTGEN
// #undef GENE_DEBUG

// 		//  if(hPt[0]*hPt[0]+hPt[1]*hPt[1] < 1.5) return PHPTMOM; // because PDF param are valid only for scale > 0.8, and it is quite safe to do this cut, as we finally cut at 2.5
// 		isub =  pypars_.msti[0];
// 		M_charm = pydat2_.pmas[0][3];
// 		int igammafrommu=3;
// 		int itarget=4;
// 		kftarget = pyjets_.k[1][itarget];
// 		xi1 = pypars_.pari[32];
// 		xi2 = pypars_.pari[33];
// 		scalePYTHIA = pypars_.pari[23];

// 		// a_ll calculation is relevant for LODIS, QCDC or PGF only

// 		igamma = 7;
// 		iparton = 8;

// 		switch(isub) {
// 			case 99:
// 				ioutq = 9;
// 				kf = pyjets_.k[1][ioutq];
// 				break;
// 			case 131:
// 			case 132:
// 				ioutq = 10;
// 				ioutq2 = 9;
// 				kf = pyjets_.k[1][ioutq];
// 				kf2 = pyjets_.k[1][ioutq2];
// 				break;
// 			case 135:
// 			case 136:
// 				ioutq = 9;
// 				ioutq2 = 10;
// 				kf = pyjets_.k[1][ioutq];
// 				kf2 = pyjets_.k[1][ioutq2];
// 				if(kf<0) {// we want the quark, not the antiquark
// 					ioutq = 10;
// 					ioutq2 = 9;
// 					kf = pyjets_.k[1][ioutq]; 
// 					kf2 = pyjets_.k[1][ioutq2]; 
// 				}
// 				break;
// 			case 11:
// 				ioutq = 9; // q coming from the photon
// 				ioutq2 = 10; // q' coming from the nucleon
// 				kf = pyjets_.k[1][ioutq];
// 				kf2 = pyjets_.k[1][ioutq2];
// 				break;
// 			case 68:
// 				ioutq = 9;
// 				ioutq2 = 10;
// 				kf = pyjets_.k[1][ioutq];
// 				kf2 = pyjets_.k[1][ioutq2];
// 				break;
// 			case 28:
// 				ioutq = 9; // parton coming from the photon
// 				ioutq2 = 10; // parton coming from the nucleon

// 				kf = pyjets_.k[1][ioutq];
// 				kf2 = pyjets_.k[1][ioutq2];
// 				break;
// 			case 13:
// 			case 53:
// 				// nothing to do here, as a_ll = -1
// 				break;
// 			case 12:
// 				ioutq = 9;
// 				ioutq2 = 10;
// 				kf = pyjets_.k[1][ioutq];
// 				kf2 = pyjets_.k[1][ioutq-2]; // we need it to calculate a_ll, see below
// 				break;
// 			default:
// 				//     assert(0);
// 				break;
// 		}
// 		if(isub<90) {
// 			kfa = pyjets_.k[1][igamma];
// 			kfb = pyjets_.k[1][iparton];
// 		}

// 		lvg = TLorentzVector(pyjets_.p[0][igammafrommu],pyjets_.p[1][igammafrommu],pyjets_.p[2][igammafrommu],pyjets_.p[3][igammafrommu]);
// 		TVector3 momgamma;
// 		momgamma = TVector3(pyjets_.p[0][igammafrommu],pyjets_.p[1][igammafrommu],pyjets_.p[2][igammafrommu]);
// 		lvt = TLorentzVector(pyjets_.p[0][itarget],pyjets_.p[1][itarget],pyjets_.p[2][itarget],pyjets_.p[3][itarget]);
// 		lvcom = lvg + lvt;
// 		comboost = lvcom.BoostVector();

// 		if( isub != 99 ) { // QCDC and PGF 
// 			lvoutq = TLorentzVector(lujets_.p[0][ioutq],lujets_.p[1][ioutq],
// 					lujets_.p[2][ioutq],lujets_.p[3][ioutq]);
// 			TVector3 momoutq;
// 			momoutq = TVector3(lujets_.p[0][ioutq],lujets_.p[1][ioutq],lujets_.p[2][ioutq]);     

// 			lvoutq2 = TLorentzVector(lujets_.p[0][ioutq2],lujets_.p[1][ioutq2],
// 					lujets_.p[2][ioutq2],lujets_.p[3][ioutq2]);
// 			TVector3 momoutq2;
// 			momoutq2 = TVector3(lujets_.p[0][ioutq2],lujets_.p[1][ioutq2],lujets_.p[2][ioutq2]); 
// 			lvingamma = TLorentzVector(lujets_.p[0][igamma],lujets_.p[1][igamma],
// 					lujets_.p[2][igamma],lujets_.p[3][igamma]);
// 			lvinparton = TLorentzVector(lujets_.p[0][iparton],lujets_.p[1][iparton],
// 					lujets_.p[2][iparton],lujets_.p[3][iparton]);

// 			TVector3 transmomoutq;
// 			TVector3 transmomoutq2;
// 			TVector3 gammadir;
// 			gammadir = momgamma.Unit();
// 			transmomoutq = momoutq - momoutq.Dot(gammadir)*gammadir;
// 			transmomoutq2 = momoutq2 - momoutq2.Dot(gammadir)*gammadir;
// 			qPt1 = sqrt(transmomoutq.Dot(transmomoutq));
// 			qPt2 = sqrt(transmomoutq2.Dot(transmomoutq2));
// 			qTheta1 = acos((gammadir.Dot(momoutq))/sqrt(momoutq.Dot(momoutq)));
// 			qTheta2 = acos((gammadir.Dot(momoutq2))/sqrt(momoutq2.Dot(momoutq2)));
// 		}



// 		// #define COMS
// #ifdef COMS
// 		// this is used for debugging purpose.
// 		// going to center of mass of gamma-p

// 		double bx=-comboost.Px();
// 		double by=-comboost.Py();
// 		double bz=-comboost.Pz();
// 		double zero=0;
// 		int izero=0;

// 		pyrobo_(&izero, &izero, &zero, &zero, &bx, &by, &bz);

// 		TLorentzVector lvg2(pyjets_.p[0][igammafrommu],pyjets_.p[1][igammafrommu],pyjets_.p[2][igammafrommu],pyjets_.p[3][igammafrommu]);
// 		double gphi = -lvg2.Phi();
// 		pyrobo_(&izero, &izero, &zero, &gphi, &zero, &zero, &zero);

// 		TLorentzVector lvg3(pyjets_.p[0][igammafrommu],pyjets_.p[1][igammafrommu],pyjets_.p[2][igammafrommu],pyjets_.p[3][igammafrommu]);
// 		double gtheta = -lvg3.Theta();
// 		pyrobo_(&izero, &izero, &gtheta, &zero, &zero, &zero, &zero);

// 		cout<<"COMS ---------"<<endl;
// 		int un=1;  
// 		pylist_(&un);
// #endif // COMS

// #endif // PHASTinPYTHIA

// #ifdef PHASTinLEPTO // ----------------------------------------
// #define INEXTGEN
// 		//  if(hPt[0]*hPt[0]+hPt[1]*hPt[1] < 0.81) return PHPTMOM;  // necessary to have scale > 0.8 for the use of GRV98 and GRSV2000
// 		// subprocess index (converting to pythia's standards)

// 		isub = leptou_.lst[23];
// 		kf = leptou_.lst[24];
// 		if(kf<0) kf = -kf;
// 		xp = leptou_.parl[27];
// 		zq = leptou_.parl[28];
// 		phi = leptou_.parl[29];
// 		M_charm = ludat2_.pmas[0][3];

// 		igamma = 2;
// 		iparton = 5;
// 		ioutq = leptou_.lst[25] - 1;
// 		ioutq2 = -1;
// 		switch(isub) {
// 			case 1:
// 				isub = 99;
// 				break;
// 			case 2:
// 				isub = 131;
// 				ioutq2 = ioutq+1;
// 				break;
// 			case 3:
// 				isub = 135;
// 				ioutq2 = ioutq+1;
// 				break;    
// 			default:
// 				cerr<<"undefined process id : "<<isub<<endl;
// 		} 

// 		int igammafrommu=igamma;
// 		int itarget=1;

// 		lvg = TLorentzVector(lujets_.p[0][igammafrommu],lujets_.p[1][igammafrommu],
// 				lujets_.p[2][igammafrommu],lujets_.p[3][igammafrommu]);
// 		lvingamma = lvg;
// 		lvt = TLorentzVector(lujets_.p[0][itarget],lujets_.p[1][itarget],
// 				lujets_.p[2][itarget],lujets_.p[3][itarget]);
// 		lvcom = lvg + lvt;
// 		comboost = lvcom.BoostVector();

// 		if( isub != 99 ) { 
// 			lvoutq = TLorentzVector(lujets_.p[0][ioutq],lujets_.p[1][ioutq],
// 					lujets_.p[2][ioutq],lujets_.p[3][ioutq]);
// 			lvoutq2 = TLorentzVector(lujets_.p[0][ioutq2],lujets_.p[1][ioutq2],
// 					lujets_.p[2][ioutq2],lujets_.p[3][ioutq2]);
// 		}    

// 		// #define COMS
// #ifdef COMS
// 		// this is used for debugging purpose.
// 		// going to center of mass of gamma-p


// 		//   TLorentzVector lvg(lujets_.p[0][2],lujets_.p[1][2],lujets_.p[2][2],lujets_.p[3][2]);
// 		//   TLorentzVector lvt(lujets_.p[0][1],lujets_.p[1][1],lujets_.p[2][1],lujets_.p[3][1]);


// 		//   TLorentzVector lvcom = lvg + lvt;
// 		//   TVector3 comboost = lvcom.BoostVector();
// 		float bx=-comboost.Px();
// 		float by=-comboost.Py();
// 		float bz=-comboost.Pz();
// 		float zero = 0;
// 		lurobo_(&zero, &zero, &bx, &by, &bz);
// 		TLorentzVector lvg2(lujets_.p[0][2],lujets_.p[1][2],lujets_.p[2][2],lujets_.p[3][2]);
// 		float gphi = -lvg2.Phi();
// 		lurobo_(&zero, &gphi, &zero, &zero, &zero);
// 		TLorentzVector lvg3(lujets_.p[0][2],lujets_.p[1][2],lujets_.p[2][2],lujets_.p[3][2]);
// 		float gtheta = -lvg3.Theta();
// 		lurobo_(&gtheta, &zero, &zero, &zero, &zero);

// 		int un=1;
// 		cout<<"COMS ------------------------------------------------"<<endl;
// 		lulist_(&un);
// #endif // COMS

// #endif // PHASTinLEPTO

// #ifdef PHASTextgenLEPTO
// #define INEXTGEN
// 		//  if(hPt[0]*hPt[0]+hPt[1]*hPt[1] < 0.81) return PHPTMOM;  // necessary to have scale > 0.8 for the use of GRV98 and GRSV2000
// 		// subprocess index (converting to pythia's standards)

// 		if(e.IsMC()) {
// 			NLUDATA ld;
// 			if(e.MCgen(ld)) {
// 				isub = ld.lst[23];
// 				kf = ld.lst[24];
// 				// if(kf<0) kf = -kf; 
// 				xp = ld.parl[27];
// 				zq = ld.parl[28];
// 				phi = ld.parl[29];

// 				igamma = 2;
// 				imu = 0;
// 				imup = 3;
// 				iparton = 5;
// 				ioutq = ld.lst[25] - 1;
// 				ioutq2 = -1;
// 				switch(isub) {
// 					case 1:
// 						isub = 99;
// 						break;
// 					case 2:
// 						isub = 131;
// 						ioutq2 = ioutq+1;
// 						break;
// 					case 3:
// 						isub = 135;
// 						if(ld.lst[7]==1) ioutq2 = ioutq+2; //PSOFF
// 						else ioutq2 = ioutq+1;
// 						break;    
// 					default:
// 						cerr<<"undefined process id : "<<isub<<endl;
// 				}

// 				int igammafrommu=igamma;
// 				int itarget=1;
// 				vector<LUJET> lujets;
// 				int nparticles = ioutq2+1;
// 				if(e.MCgen(nparticles,lujets)) {
// 					lvmu_g = TLorentzVector(lujets[imu].p[0],lujets[imu].p[1],
// 							lujets[imu].p[2],lujets[imu].p[3]);
// 					lvmup_g = TLorentzVector(lujets[imup].p[0],lujets[imup].p[1],
// 							lujets[imup].p[2],lujets[imup].p[3]);
// 					lvg = TLorentzVector(lujets[igammafrommu].p[0],lujets[igammafrommu].p[1],
// 							lujets[igammafrommu].p[2],lujets[igammafrommu].p[3]);
// 					lvingamma = lvg;
// 					lvt = TLorentzVector(lujets[itarget].p[0],lujets[itarget].p[1],
// 							lujets[itarget].p[2],lujets[itarget].p[3]);
// 					lvcom = lvg + lvt;
// 					comboost = lvcom.BoostVector();

// 					if( isub != 99 ) { 
// 						// Energy and momentum returned by Lepto do not give proper
// 						// mass. Calculate energy from momentum and mass.
// 						double energy1 = sqrt(lujets[ioutq].p[0]*lujets[ioutq].p[0] + lujets[ioutq].p[1]*lujets[ioutq].p[1] +
// 								lujets[ioutq].p[2]*lujets[ioutq].p[2] + lujets[ioutq].p[4]*lujets[ioutq].p[4]);
// 						// double energy1 = lujets[ioutq].p[3]; 
// 						lvoutq = TLorentzVector(lujets[ioutq].p[0],lujets[ioutq].p[1],
// 								lujets[ioutq].p[2],energy1);
// 						kf = lujets[ioutq].k[1];
// 						double energy2 = sqrt(lujets[ioutq2].p[0]*lujets[ioutq2].p[0] + lujets[ioutq2].p[1]*lujets[ioutq2].p[1] +
// 								lujets[ioutq2].p[2]*lujets[ioutq2].p[2] + lujets[ioutq2].p[4]*lujets[ioutq2].p[4]);
// 						// double energy2 = lujets[ioutq2].p[3]; 
// 						lvoutq2 = TLorentzVector(lujets[ioutq2].p[0],lujets[ioutq2].p[1],
// 								lujets[ioutq2].p[2],energy2);
// 						kf2 = lujets[ioutq2].k[1];

// #ifdef GENE_DEBUG // Set it in HptEventOld.h
// 						event->ggM = lvg.M();
// 						event->gqM1 = lvoutq.M();
// 						event->gqM2 = lvoutq2.M();
// 						event->ggE = lvg.E();
// 						event->gqE1 = lvoutq.E();
// 						event->gqE2 = lvoutq2.E();
// #endif //GENE_DEBUG
// 					}

// #ifdef GENE_DUMP
// 					cout<<"\n############"<<endl;
// 					cout<<"R: "<<eRun<<" S: "<<eSpill<<" E: "<<eEvinspill<<endl;
// 					cout<<"Process: "<<isub<<endl;
// 					cout<<"ID:ParticleCode:Parent:Doughter:P:E:M:px:py:pz"<<endl;
// 					for(int i=0; i<nparticles; i++) {
// 						cout<<i+1<<'\t'<<lujets[i].k[1]<<'\t'
// 							<<lujets[i].k[2]<<'\t'<<lujets[i].k[3]<<'\t'
// 							<<sqrt(lujets[i].p[0]*lujets[i].p[0] + lujets[i].p[1]*lujets[i].p[1] + lujets[i].p[2]*lujets[i].p[2])<<'\t'
// 							<<lujets[i].p[3]<<'\t'<<lujets[i].p[4]<<'\t'
// 							<<lujets[i].p[0]<<'\t'<<lujets[i].p[1]<<'\t'<<lujets[i].p[2]<<'\t'
// 							<<endl;
// 					}
// #endif //GENE_DUMP
// 				}
// 			}
// 		}
// #endif // #ifdef PHASTextgenLEPTO


// #ifdef INEXTGEN //---------------------------------------------

// 		// calculation of a_ll for lepto or pythia. 

// 		float xptmp=-1;
// 		double xbjtmp=-1;
// 		double nutmp=-1;
// 		double q2tmp=-1;
// 		float zqtmp=-1;
// 		float zqtmp2=-1;
// 		float phitmp=-1;
// 		float phitmp2=-1;

// #ifdef PHASTextgenLEPTO
// 		nutmp = lvmu_g.E() - lvmup_g.E();
// #else
// 		nutmp = pq/M_p;
// #endif
// 		q2tmp = -lvg.M2();
// 		xbjtmp = q2tmp / (2.0 * M_p * nutmp);

// 		if(isub != 99) {
// 			TLorentzVector lvhat = lvoutq + lvoutq2;
// 			TLorentzVector lvthat = lvingamma - lvoutq;
// 			TLorentzVector lvuhat = lvingamma - lvoutq2;

// 			shat = lvhat.M2();
// 			that = lvthat.M2();
// 			uhat = lvuhat.M2();
// 		}
// 		if(isub!=99) {
// #ifdef GENE_DEBUG // Set it in HptEventOld.h
// 			gKsi_dbg = (shat+kQ2)/(2*pq);
// #endif //GENE_DEBUG
// 			ksi = (shat+q2tmp)/(2*M_p*nutmp);

// #ifdef GENE_DEBUG // Set it in HptEventOld.h
// 			gXp_dbg = kXbj/ksi; 
// #endif //GENE_DEBUG
// 			xptmp = xbjtmp/ksi;

// 			// calculation of zq
// 			// Standard method based on definition in POLDIS paper
// 			// zqtmp = lvt.Dot(lvoutq)/lvt.Dot(lvg);
// 			// Naive ;) method - works as fine as standard. No problem with not well defined 4mom of target.
// 			zqtmp = lvoutq.E() / nutmp;
// 			zqtmp2 = lvoutq2.E() / nutmp;
// #ifdef GENE_DEBUG // Set it in HptEventOld.h
// 			// Calculate z_q from uhat. This is based on POLDIS appendix.
// 			gZq_dbg = (lvoutq.M2() - uhat) * (xptmp / q2tmp);
// 			// Calculate z_qhat
// 			gZq2_dbg = (lvoutq2.M2() - that) * (xptmp / q2tmp);
// #endif //GENE_DEBUG

// 			// calculation of phi 
// #ifdef GENE_DEBUG // Set it in HptEventOld.h
// 			TVector3 v3perp2gammaq = v3q.Cross( lvoutq.Vect() );
// 			TVector3 v3perp2gammaq2 = v3q.Cross( lvoutq2.Vect() );
// 			v3perp2gammaq = v3perp2gammaq.Unit();
// 			v3perp2gammaq2 = v3perp2gammaq2.Unit();
// 			if( v3q.Dot( v3perp2qdiff.Cross(v3perp2gammaq) ) >= 0)
// 				gPhiq_dbg = acos(v3perp2gammaq * v3perp2qdiff);
// 			else 
// 				gPhiq_dbg = -acos(v3perp2gammaq * v3perp2qdiff) + 2*TMath::Pi();     

// 			if( v3q.Dot( v3perp2qdiff.Cross(v3perp2gammaq2) ) >= 0)
// 				gPhiq2_dbg = acos(v3perp2gammaq2 * v3perp2qdiff);
// 			else 
// 				gPhiq2_dbg = -acos(v3perp2gammaq2 * v3perp2qdiff) + 2*TMath::Pi();     
// #endif //GENE_DEBUG

// 			// New calculation of phi by Konrad
// 			// Phiq is an azimuthal angle of outgoing quark defined wrt lepton
// 			// scattering plane in gamma-proton CMS.  Actually for formulas derived
// 			// by Krzysztof it is difference of azimuthal angles of outgoing quark
// 			// and lepton in Breit system.  gP CMS is a good approximation of
// 			// Breit. X axis can be choosen so that lepton azimuth angle is zero so
// 			// those two definitions are compatible.  In gP CMS Z axis is defined
// 			// as direction of the photon. Then XY plane is a plane that gamma is
// 			// normal to.  In three-dimensional polar coordinate systems, including
// 			// cylindrical coordinates and spherical coordinates, the azimuth of a
// 			// point is the angle between the positive x-axis and the projection of
// 			// the vector onto the xy-plane (the component of the vector in the
// 			// xy-plane).
// 			//
// 			// Prepare the vectors.
// 			TLorentzVector lvg_b(lvg);
// 			TLorentzVector lvmu_b(lvmu);
// 			TLorentzVector lvmup_b(lvmup);
// 			TLorentzVector lvoutq_b(lvoutq);
// 			TLorentzVector lvoutq2_b(lvoutq2);

// 			// Boost to CMS of P and g
// 			lvg_b.Boost(-comboost);
// 			lvmu_b.Boost(-comboost);
// 			lvmup_b.Boost(-comboost);
// 			lvoutq_b.Boost(-comboost);
// 			lvoutq2_b.Boost(-comboost);

// 			// Get 3 vectors in CMS of P and g
// 			TVector3 v3mu_b = lvmu_b.Vect();
// 			TVector3 v3mup_b = lvmup_b.Vect();
// 			TVector3 v3outq_b = lvoutq_b.Vect();
// 			TVector3 v3outq2_b = lvoutq2_b.Vect();
// 			// Gamma defines Z direction
// 			// We will project everything to a plane 0XY
// 			TVector3 v3norm = lvg_b.Vect().Unit();
// 			// Project quark and leptons to 0XY plane -> P(q), P(mu), P(mu')
// 			// P(A) = Norm x (A x Norm)
// 			v3outq_b = v3outq_b.Cross(v3norm);
// 			v3outq_b = v3norm.Cross(v3outq_b);
// 			v3outq2_b = v3outq2_b.Cross(v3norm);
// 			v3outq2_b = v3norm.Cross(v3outq2_b);
// 			v3mu_b = v3mu_b.Cross(v3norm);
// 			v3mu_b = v3norm.Cross(v3mu_b);
// 			v3mup_b = v3mup_b.Cross(v3norm);
// 			v3mup_b = v3norm.Cross(v3mup_b);
// 			// Calculate angle between P(q) and P(mu)
// 			// We want to get the azimuthal angle (0, 2Pi) calculated from P(mu).
// 			// X axis is defined by P(mu) direction.  By checking if product of
// 			// P(mu) x P(q) points in Z direction we check in which half of XY
// 			// plane P(q) lays. Y>0 or Y<0.  With this knowladge we can use
// 			// T3Vector::Angle to calculate the angle.
// 			TVector3 rot_dir = v3mu_b.Cross(v3outq_b).Unit();
// 			if(v3norm.Angle(rot_dir) < 0.5*TMath::Pi()) phitmp = v3mu_b.Angle(v3outq_b);
// 			else phitmp = 2*TMath::Pi() - v3mu_b.Angle(v3outq_b);

// 			rot_dir = v3mu_b.Cross(v3outq2_b).Unit();
// 			if(v3norm.Angle(rot_dir) < 0.5*TMath::Pi()) phitmp2 = v3mu_b.Angle(v3outq2_b);
// 			else phitmp2 = 2*TMath::Pi() - v3mu_b.Angle(v3outq2_b);

// #ifdef GENE_DEBUG // Set it in HptEventOld.h
// 			event->gPhi_lepton = v3mu_b.Angle(v3mup_b); // Store angle between P(mu) and P(mu'). It should be 0.
// #endif //GENE_DEBUG
// 		} else if(isub==99) {
// 			ksi = 1;
// 			xi1=1;
// 			xi2=event->kXbj;
// 			xptmp = 1;
// 			zqtmp = 1; 
// 			zqtmp2 = 0; 
// 			phitmp = 0;
// 			phitmp2 = 0;
// 		}

// #ifdef GENE_DUMP2
// 		cout<<"****"<<endl;
// 		cout<<"gamma:kf, P, M, E, px, py, pz"<<endl;
// 		cout<<22<<'\t'<<lvg.P()<<'\t'<<lvg.M()<<'\t'<<lvg.E()<<'\t'<<lvg.Px()<<'\t'<<lvg.Py()<<'\t'<<lvg.Pz()<<endl;
// 		cout<<"First quark:kf, P, M, E, px, py, pz"<<endl;
// 		cout<<kf<<'\t'<<lvoutq.P()<<'\t'<<lvoutq.M()<<'\t'<<lvoutq.E()<<'\t'<<lvoutq.Px()<<'\t'<<lvoutq.Py()<<'\t'<<lvoutq.Pz()<<endl;
// 		cout<<"Second quark:kf, P, M, E"<<endl;
// 		cout<<kf2<<'\t'<<lvoutq2.P()<<'\t'<<lvoutq2.M()<<'\t'<<lvoutq2.E()<<'\t'<<lvoutq2.Px()<<'\t'<<lvoutq2.Py()<<'\t'<<lvoutq2.Pz()<<endl;
// 		cout<<"shat, that, uhat, shat+that+uhat, 2*m_Q^2 - Q^2"<<endl;
// 		cout<<shat<<'\t'<<that<<'\t'<<uhat<<'\t'<<shat+that+uhat<<'\t'<<2*lvoutq.M2()-q2tmp<<endl;
// 		cout<<"xp, zq, phi (Lepto gen)"<<endl;
// 		cout<<xp<<'\t'<<zq<<'\t'<<phi<<endl;
// 		cout<<"xp, zq1, phi1, zq2, phi2, zq1+zq2, phi_lepton (calc)"<<endl;
// 		cout<<xptmp<<'\t'<<zqtmp<<'\t'<<phitmp<<'\t'<<zqtmp2<<'\t'<<phitmp2<<'\t'<<zqtmp+zqtmp2<<'\t'<<gPhi_lepton<<endl;
// #endif //GENE_DUMP2

// 		// otherwise calculated variables are taken for a_LL calculation, and stored
// 		// in the tree at the usual place 
// 		// ie LEPTO in non debug mode, or PYTHIA
// 		event->gXp_gen = xp;
// 		event->gZq_gen = zq; 
// 		event->gPhiq_gen = phi;
// 		xp = xptmp;
// 		zq = zqtmp;
// 		zq2 = zqtmp2;
// 		phi = phitmp;
// 		phi2 = phitmp2;

// #ifndef NO_FORTRAN
// 		switch(isub) {
// 			case 99:
// 				all = allq_(&(event->kY));	            
// 				event->gAll_gen = all;
// 				break;
// 			case 131:
// 			case 132:
// 				//    all = allqg_(&kY, &xp, &zq, &phi);
// 				//all = allqgkurek_(&(event->kY), &xp, &zq, &(event->kQ2), &(event->kXbj),&phi);	      
// 				//event->gAll_gen = allqgkurek_(&(event->kY), &(event->gXp_gen), &(event->gZq_gen),&(event->kQ2), &(event->kXbj),&(event->gPhiq_gen));	      
// 				break;
// 			case 135:
// 			case 136:
// 				if(abs(kf) == 4) {
// 					//     all = allqqhf_(&kY, &kQ2, &xp, &zq, &phi, &M_charm);
// 					//all = allqqhfkurek_(&(event->kY), &xp, &zq, &(event->kQ2), &(event->kXbj), &M_charm, &phi);
// 					event->gAll_gen = allqqhfkurek_(&(event->kY), &(event->gXp_gen), &(event->gZq_gen),
// 							&(event->kQ2), &(event->kXbj), &M_charm, &(event->gPhiq_gen));	      
// 				} else {
// 					//    all = allqq_(&kY, &xp, &zq, &phi); 
// 					all = allqqkurek_(&(event->kY), &xp, &zq, &(event->kQ2), &(event->kXbj), &phi);
// 					event->gAll_gen = allqqkurek_(&(event->kY), &(event->gXp_gen), &(event->gZq_gen),
// 							&(event->kQ2), &(event->kXbj), &(event->gPhiq_gen));	      
// 				}
// 				break;
// 			case 68:
// 				all = all68_(&shat, &that, &uhat);
// 				break;
// 			case 28:
// 				all = all28_(&shat, &that, &uhat);
// 				break;
// 			case 11:
// 				if(kf == kf2) all = all11qq_(&shat, &that, &uhat);
// 				if(kf == -kf2) all = all12qqbar_(&shat, &that, &uhat); 
// 				if(abs(kf)!=abs(kf2)) all = all11qqprim_(&shat, &that, &uhat);
// 				break;
// 			case 13:
// 			case 53:
// 				all = -1;
// 				break;
// 			case 12:
// 				all = -1;
// 				break;
// 			default:
// 				break;
// 		}

// 		// filling uservar array. first 3 variables are here just in case of...
// 		// if mods are introduced here, modify also the code for the reading mode
// 		// below

// 		//cout<<isub<<"\t"<<xi1<<"\t"<<xi2<<"\t"<<kXbj<<"\t"<<xp<<endl;

// #ifndef PHASTextgenLEPTO
// 		uservar_.uv[3]  = event->kQ2;
// 		uservar_.uv[4]  = event->kXbj;
// 		uservar_.uv[5]  = event->kY;
// 		uservar_.uv[6]  = shat;
// 		uservar_.uv[7]  = uhat;
// 		uservar_.uv[8]  = that;
// 		uservar_.uv[9]  = xp;
// 		uservar_.uv[10] = zq;
// 		uservar_.uv[11] = phi;
// 		uservar_.uv[12] = phi2;
// 		uservar_.uv[13] = qTheta1;
// 		uservar_.uv[14] = qTheta2;
// 		uservar_.uv[15] = all;
// 		uservar_.uv[19] = e.RunNum();

// 		// -- calculation of dp/p for parton in nucleon, resolved process

// 		// if(nhadrons>1)    double scale = 1.*(hPt[0]*hPt[0]+hPt[1]*hPt[1]); // choice of scale !
// 		double scale = 0.;
// 		for(int i = 0;i<min(hN,2);i++) scale = scale + event->hPt[i]*event->hPt[i];
// 		//    double scale = qPt1*qPt1+qPt2*qPt2;
// 		scale = max(scale,event->kQ2);
// 		// double scale = 10; // just to ensure > 0.8; for tests without any cuts
// 		cout << "scale is: " << scale << endl;
// 		cout << "Pt0 is: " << event->hPt[0] << endl;
// 		cout << "Pt1 is: " << event->hPt[1] << endl;
// 		//  if(scale<0.81) {
// 		//    cout << "Aaaargh!! scale is too low=> i reject this event" << endl;
// 		//    return PHPTMOM;
// 		// }
// 		double P2 = event->kQ2;
// 		//  if(P2>0.5*scale) {
// 		//    cout << "have to ignore this event, P2>0.5*scale" << endl;
// 		//return PHPTMOM2;
// 		// provisoire!!!
// 		//  }

// 		if(isub<90 && P2<0.5*scale && scale>0.81) {

// 			double uv , dv , us , ds , ss , gl, u , d , ub , db , st , gl2 , g1p , g1n;
// 			int iset = 1;
// 			int iset2 = 3;

// 			// cout << "Pt0 is: " << hPt[0] << endl;
// 			// cout << "Pt1 is: " << hPt[1] << endl;
// 			grv98pa_(&iset, &xi2, &scale, &uv, &dv, &us, &ds, &ss, &gl);
// 			parpol_(&iset2, &xi2, &scale, &u, &d, &ub, &db, &st, &gl2, &g1p, &g1n);
// 			switch(kf2) {
// 				case -1: // parton is dbar
// 					if(kftarget==2212) delpovpnucleon = db/ds;	// target is proton
// 					if(kftarget==2112) delpovpnucleon = ub/us;        // target is neutron
// 					break;
// 				case -2: // parton is ubar
// 					if(kftarget==2212) delpovpnucleon = ub/us;
// 					if(kftarget==2112) delpovpnucleon = db/ds;            
// 					break;
// 				case 2: // parton is u
// 					if(kftarget==2212) delpovpnucleon = u/(uv+us);
// 					if(kftarget==2112) delpovpnucleon = d/(dv+ds);            
// 					break;
// 				case 1: // parton is d
// 					if(kftarget==2212) delpovpnucleon = d/(dv+ds);
// 					if(kftarget==2112) delpovpnucleon = u/(uv+us);   
// 					break;
// 				case 3:
// 				case -3: // parton is s or sbar
// 					delpovpnucleon = st/ss;	            
// 					break;
// 				case 21:
// 					delpovpnucleon = 0; //  in this case, deltaG/G is involved!
// 				default:
// 					break;
// 			}

// 			cout << "delpovpnucleon = " <<  delpovpnucleon << endl;
// 			// -- end of calculation of dp/p for p in nucleon

// 			// -- now the same but for p in photon
// 			// -- case 1: photon is real
// 			double uphmax, dphmax, sphmax, cphmax, glphmax;
// 			double uphmin, dphmin, sphmin, cphmin, glphmin;
// 			int imin = 2; // min scenario
// 			int imax = 1; // max scenario
// 			parphot_(&imin,&xi1,&scale,&uphmin, &dphmin, &sphmin, &cphmin, &glphmin);
// 			parphot_(&imax,&xi1,&scale,&uphmax, &dphmax, &sphmax, &cphmax, &glphmax);
// 			// -- case 2: photon is virtual
// 			// -- part 1: polarized PDF
// 			double uvphmax, dvphmax, svphmax, glvphmax;
// 			double uvphmin, dvphmin, svphmin, glvphmin;
// 			//  double P2 = kQ2;
// 			////  int iminv = 2; // min scenario
// 			////  int imaxv = 1; // max scenario
// 			//  if(P2>0.5*scale) {
// 			//    cout << "have to ignore this event, P2>0.5*scale" << endl;
// 			//return PHPTMOM2;
// 			// provisoire!!!
// 			//  }
// 			parvpol_(&imax,&xi1,&scale,&P2, &uvphmax, &dvphmax, &svphmax, &glvphmax);
// 			parvpol_(&imin,&xi1,&scale,&P2, &uvphmin, &dvphmin, &svphmin, &glvphmin);


// 			// -- part 2: unpolarized PDF 
// 			double uphv, dphv, sphv, gphv;
// 			grsglo_(&xi1, &scale, &P2, &uphv, &dphv, &sphv, &gphv);

// 			//  double s, alphaprim, beta, a,b, AA, BB, CC, DD, EE, Eprim, sf;
// 			//  s = (log((log(scale/(0.232*0.232)))/log(0.25/(0.232*0.232))))/log(2.718);
// 			//  double unpolpphoton;

// 			switch(kf) {
// 				case 2:
// 				case -2:
// 					//    alphaprim = 1.717; beta = 0.641; a = 0.5-0.176*s;
// 					//   b = 15.00-5.687*sqrt(s) - 0.552*s*s;
// 					//  AA = 0.235+0.046*sqrt(s);
// 					//  BB = 0.082-0.051*s + 0.168*s*s;
// 					//  CC = 0.459*s; DD = 0.354 - 0.061*s;
// 					//  EE = 4.899 + 1.678*s; Eprim = 2.046 + 1.389*s;

// 					//    unpolpphoton = unpolphotudg_(&xi1,&AA,&BB,&CC,&DD,&EE,&Eprim,&a,&b,&s,&alphaprim,&beta);

// 					delpovpgammin = uvphmin / uphv;
// 					delpovpgammax = uvphmax / uphv;

// 					break;

// 				case 1:
// 				case -1:
// 					//   alphaprim = 1.549; beta = 0.782; a = 0.496+0.026*s;
// 					//  b = 0.685-0.580*sqrt(s) + 0.608*s*s;
// 					// AA = 0.233+0.302*s;
// 					// BB = -0.818*s + 0.198*s*s;
// 					// CC = 0.114 + 0.154*s; DD = 0.405 - 0.195*s + 0.046*s*s;
// 					// EE = 4.807 + 1.226*s; Eprim = 2.166 + 0.664*s;

// 					//  unpolpphoton = unpolphotudg_(&xi1,&AA,&BB,&CC,&DD,&EE,&Eprim,&a,&b,&s,&alphaprim,&beta);

// 					delpovpgammin = dvphmin / dphv;
// 					delpovpgammax = dvphmax / dphv;

// 					break;
// 				case 3:
// 				case -3:
// 					//   sf = 0; alphaprim = 1.609; beta = 0.962; a = 0.470 - 0.099*s*s;
// 					//  b = 3.246;
// 					// AA = 0.121 - 0.068*sqrt(s);
// 					// BB = -0.090 + 0.074*s;
// 					// CC = 0.062 + 0.034*s; DD = 0.226*s - 0.060*s*s;
// 					// EE = 4.288 + 1.707*s; Eprim = 2.122 + 0.656*s;

// 					//  unpolpphoton = unpolphotscb_(&xi1,&AA,&BB,&CC,&DD,&EE,&Eprim,&a,&b,&s,&alphaprim,&beta, &sf);

// 					delpovpgammin = svphmin / sphv;
// 					delpovpgammax = svphmax / sphv;

// 					break;
// 				case 21:
// 					//   alphaprim = 0.676; beta = 1.089; a = 0.462 - 0.524*sqrt(s);
// 					//  b = 5.451 - 0.804*s*s;
// 					// AA = 0.535 - 0.504*sqrt(s) + 0.288*s*s;
// 					//  BB = 0.364 - 0.520*s;
// 					// CC = -0.323 + 0.115*s*s; DD = 0.233 + 0.790*s - 0.139*s*s;
// 					// EE = 0.893 + 1.968*s; Eprim = 3.432 + 0.392*s;

// 					//  unpolpphoton = unpolphotudg_(&xi1,&AA,&BB,&CC,&DD,&EE,&Eprim,&a,&b,&s,&alphaprim,&beta);

// 					delpovpgammin = glvphmin / gphv;
// 					delpovpgammax = glvphmax / gphv;

// 					break;

// 				default:
// 					break;
// 			}

// 			cout <<"delpovpgammin = " << delpovpgammin << endl; 
// 			cout <<"delpovpgammax = " << delpovpgammax << endl; 

// /*
// 			uservar_.uv[16] = delpovpnucleon;   
// 			uservar_.uv[17] = delpovpgammin;   
// 			uservar_.uv[18] = delpovpgammax;  
// 			uservar_.uv[19]  = kf2; 
// */
// 			//  uservar_.uv[21]  = that;
// 			//  uservar_.uv[22]  = uhat;

// 			//  cout << "qPt1 = " << qPt1 << endl;
// 		}
// #endif // #ifndef PHASTextgenLEPTO
// #endif // #ifdef NO_FORTRAN

// #else // #ifdef INEXTGEN

// 		// reading mode. In case of MC data, contents of the uservar array are 
// 		// filled in the tree

// 		int impmcgen = 0;
// 		if(e.IsMC()) {
// 			impmcgen = 1;
// 			NLUDATA ld;
// 			if(e.MCgen(ld)) {
// 				impmcgen = 2;

// 				// kQ2  
// 				// kXbj 
// 				// kY 
// 				// ... are calculated from the reconstructed tracks.
// 				isub             = ld.lst[23]; //first Order QCD process
// 				kf               = ld.lst[24]; //struck quark
// 				shat             = ld.uservar[6];
// 				uhat             = ld.uservar[7];
// 				that             = ld.uservar[8];
// 				event->gXp_gen   = ld.parl[27];
// 				xp               = ld.uservar[9];
// 				event->gZq_gen   = ld.parl[28];
// 				zq               = ld.uservar[10];
// 				event->gPhiq_gen = ld.parl[29];
// 				phi              = ld.uservar[11];
// 				phi2             = ld.uservar[12];
// 				qTheta1          = ld.uservar[13];
// 				qTheta2          = ld.uservar[14];
// 				all              = ld.uservar[15];
// 				event->eEvNr     = e.RunNum()*10000 + ld.uservar[19];
// 			}
// 		} else {
// 			event->eEvNr = e.UniqueEvNum();
// 		}

// #endif
// 		if(isub==28 && kf>5) event->gIsub = isub+1; // now gq->gq is process 29
// 		else  event->gIsub = isub;

// 		event->gKf = kf;
// 		if(isub!=12) event->gKf2 = kf2;
// 		event->gShat = shat;
// 		event->gThat = that;
// 		event->gUhat = uhat;
// 		event->gKsi=ksi;
// 		event->gXp = xp;
// 		event->gZq = zq;
// 		event->gZq2 = zq2;
// 		event->gPhiq = phi;
// 		event->gPhiq2 = phi2;
// 		event->gAll = all;
// 		// gAllK = allK;


// 		// Those varaibles are only calculated so far for PYTHIA
// 		// TODO check if we could calculate them for lepto and if they are needed.
// #ifdef PHASTinPYTHIA
// 		// event->gXi1 = xi1; 
// 		// event->gScale = scalePYTHIA; 
// 		// event->gqPt1 = qPt1; 
// 		// event->gqPt2 = qPt2; 
// 		// event->gqTheta1 = qTheta1; 
// 		// event->gqTheta2 = qTheta2; 
// 		// event->gXi2 = xi2; 
// 		// event->gDelpovpnucleon = delpovpnucleon; 
// 		// event->gDelpovpgammin = delpovpgammin; 
// 		// event->gDelpovpgammax = delpovpgammax; 
// 		if(event->gIsub!=99){
// 			//  cout << "xi1 = " << gXi1 << endl;
// 			//  cout << "xi2 = " << gXi2 << endl;
// 			//    cout << "kQ2 = " << kQ2 << endl;
// 			//   cout << "sum(P_T^2) = " << hPt[0]*hPt[0]+hPt[1]*hPt[1] << endl;
// 			//    cout <<"gIsub = " << gIsub << endl;
// 			//  cout << "gKf = " << gKf << endl;
// 			//  cout << "gKf2 = " << gKf2 << endl;
// 			//  cout << "a_ll= " <<gAll<<endl;
// 			//  cout << "delpovpnucleon = " << gDelpovpnucleon << endl;
// 			//    cout << "------------------" << endl;
// 		}
// #endif //PHASTinPYTHIA
// #ifdef PHASTextgenLEPTO
// 		event->eEvNr = e.UniqueEvNum();
// #endif //PHASTextgenLEPTO
// 	} // END Generator block

	t->Fill(); 
	Phast::Ref().h_file = t->GetCurrentFile(); 

	e.TagToSave();
	event->Clear();
	//cout<< e.UniqueEvNum()<<endl;
	return PHPTOK;
}

