#ifndef __UTILS_H__
#define __UTILS_H__
#include "TROOT.h"
#include "PaEvent.h"

bool isVetoFired(int eTrig, Double_t Vp, Double_t Vtot, Double_t Vcalo);
int  IsCrossHodoscopes(const PaTrack& TrackMu1, int TMsk, bool Hadron);
double DepolFactor(double x, double y, double q2);
void RemoveBadMiddleTrigger(PaEvent& e, const PaTrack& mup_tr);

#endif


