// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "veto_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

/*BDX-MINI*/
double att_l_ang_IV = 120;
double att_l_ang_OV = 140;
double TrGrvIL = 0.6;
double TrGrvOL = 0.65;


double MeV2peIV[8] = { 57., 57., 66.6, 52.6, 50., 57.1, 66.6, 50. };
double MeV2peOV[8] = { 38.4615, 33.3333, 200, 40, 33.8983, 42.5532, 40, 41.6667 };

//tuning
double tunOV[8] = {1.14, 1.14, 1., 1.08, 1.2, 1.23, 1.23, 1.23};
double tunIV[8] = {1.14, 1.08, 1.14, 1.14, 1.2, 1.08, 1.08, 1.14};



double LightSpeedAng = 15. * 360 / (2 * 3.1415 * 9.);

double att_z_IV = 100 * cm;
double att_z_OV = 100 * cm;

typedef struct groupedHits {
	double Esum;
	double avgT;
	G4ThreeVector avgPos;
	G4ThreeVector avgLPos;
};

typedef struct aStep {
	bool isMatchedToGroupedHit;
	double E, T;
	G4ThreeVector pos;
	G4ThreeVector lpos;
};

//sort in DESCENDING energy
bool compareSteps(const aStep &a, const aStep &b) {
	return a.E > b.E;
}

bool isMatched(const aStep &a, const aStep &b, double dT = 1 * ns, double dX = 1 * cm) {
	bool ret = false;
	if ((fabs(a.T - b.T) <= dT) && ((a.pos - b.pos).mag() < dX)) {
		ret = true;
	}
	return ret;
}

/*aHit contains all the steps in the IV or OV for BDXmini, within a large TimeWindow(TW).
 This methods groups them together according to time and distance (dT=1 ns, dX=1cm)
 */
vector<groupedHits> sortBDXMiniVetoHits(MHit* aHit, double dT = 1 * ns, double dX = 1 * cm) {
	vector<groupedHits> ghits;

	auto times = aHit->GetTime();
	auto pos = aHit->GetPos();
	auto lpos = aHit->GetLPos();
	auto ene = aHit->GetEdep();

	vector<vector<int>> groupIDs;
	vector<int> thisIDs;

	//sort in deposited energy keeping order
	aStep step;
	vector<aStep> steps;
	for (int ii = 0; ii < pos.size(); ii++) {
		step.E = ene[ii];
		step.T = times[ii];
		step.pos = pos[ii];
		step.lpos = lpos[ii];
		step.isMatchedToGroupedHit = false;
		steps.push_back(step);
	}
	std::sort(steps.begin(), steps.end(), compareSteps);

	for (int ii = 0; ii < steps.size(); ii++) {
		if (steps[ii].isMatchedToGroupedHit == false) { //create a new group
			thisIDs.clear();
			steps[ii].isMatchedToGroupedHit = true;
			thisIDs.push_back(ii);

			for (int jj = ii + 1; jj < steps.size(); jj++) {
				if (steps[jj].isMatchedToGroupedHit == false) {
					bool match = isMatched(steps[ii], steps[jj]);
					if (match) {
						steps[jj].isMatchedToGroupedHit = true;
						thisIDs.push_back(jj);
					}
				}
			}
			groupIDs.push_back(thisIDs);
		}
	}

	for (auto group : groupIDs) {
		groupedHits ghit;
		ghit.Esum = 0;
		ghit.avgT = 0;
		ghit.avgPos.set(0, 0, 0);
		ghit.avgLPos.set(0, 0, 0);
		for (auto id : group) {
			ghit.Esum += steps[id].E;
			ghit.avgT += steps[id].T * steps[id].E;
			ghit.avgPos.setX(ghit.avgPos.x() + steps[id].E * steps[id].pos.x());
			ghit.avgPos.setY(ghit.avgPos.y() + steps[id].E * steps[id].pos.y());
			ghit.avgPos.setZ(ghit.avgPos.z() + steps[id].E * steps[id].pos.z());
			ghit.avgLPos.setX(ghit.avgLPos.x() + steps[id].E * steps[id].lpos.x());
			ghit.avgLPos.setY(ghit.avgLPos.y() + steps[id].E * steps[id].lpos.y());
			ghit.avgLPos.setZ(ghit.avgLPos.z() + steps[id].E * steps[id].lpos.z());
		}
		if (ghit.Esum > 0) {
			ghit.avgT /= ghit.Esum;
			ghit.avgPos.setX(ghit.avgPos.x() / ghit.Esum);
			ghit.avgPos.setY(ghit.avgPos.y() / ghit.Esum);
			ghit.avgPos.setZ(ghit.avgPos.z() / ghit.Esum);

			ghit.avgLPos.setX(ghit.avgLPos.x() / ghit.Esum);
			ghit.avgLPos.setY(ghit.avgLPos.y() / ghit.Esum);
			ghit.avgLPos.setZ(ghit.avgLPos.z() / ghit.Esum);

			ghits.push_back(ghit);
		}
	}

	/*	std::cout << "STEPS: " << pos.size() << std::endl;
	 int aa = 0;
	 for (auto step : steps) {
	 std::cout << aa << " : " << step.E << " " << step.T << " " << step.pos.x() << " " << step.pos.y() << " " << step.pos.z() << std::endl;
	 aa++;
	 }
	 std::cout << "GROUPS: " << groupIDs.size() << std::endl;
	 aa = 0;
	 for (auto group : groupIDs) {
	 std::cout << aa << " ";
	 for (auto id : group)
	 std::cout << id << " ";
	 std::cout << std::endl;
	 aa++;
	 }
	 std::cout << "GROUPED HITS :" << ghits.size() << std::endl;
	 aa = 0;
	 for (auto ghit : ghits) {
	 std::cout << aa << " " << ghit.Esum << " " << ghit.avgT << " " << ghit.avgPos.x() << " " << ghit.avgPos.y() << " " << ghit.avgPos.z() << std::endl;
	 aa++;
	 }*/

	return ghits;
}

map<string, double> veto_HitProcess::integrateDgt(MHit* aHit, int hitn) {
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector = identity[0].id;
	int veto_id = identity[1].id;
	int channel = identity[2].id;

	// Digitization Parameters

// initialize ADC and TDC
	double etotL = 0;
	double etotR = 0;
	double timeL = 0;
	double timeR = 0;
	int ADC1 = 0;
	int ADC2 = 0;
	int ADC3 = 0;
	int ADC4 = 0;
	int ADC5 = 0;
	int ADC6 = 0;
	int ADC7 = 0;
	int ADC8 = 0;
	int TDC1 = 4096;
	int TDC2 = 4096;
	int TDC3 = 4096;
	int TDC4 = 4096;
	int TDC5 = 4096;
	int TDC6 = 4096;
	int TDC7 = 4096;
	int TDC8 = 4096;

	double length;
	// From measurement
	// spe = 0.36pC
	// Cosmic = 84pC (~235pe)
	// Attenuation at 80cm: ~0.85 -> effective att lenght ~350cm
	//Plastic
	double paddle_surface;        // paddle surface
	double light_yield;  // number of optical photons pruced in the scintillator per MeV of deposited energy
	double att_length;               // light at tenuation length
	double veff;            // light velocity in scintillator
	//PMT
	double sensor_surface;   // area of photo sensor
	double sensor_effective_area; // considering only a fraction of the photocathod
	double sensor_qe;                     // photo sensor quantum efficiency
	double sensor_gain;         // pmt gain x electron charge in pC (2.2x10^6)x(1.6x10^-7) -> ~0.36pC or 1 to have pe
	double light_coll; // ratio of photo_sensor area over paddle section ~ light collection efficiency
	double light_guide_att;
	double tL;
	double tR;
	double peL = 0.;
	double peR;
	double etot_g4 = 0.;

	// Proposal
	// sector: run over channels, from 0 to N
	// Oveto == 5
	// channel: run over the position (1=T 2=B 3=R 4=L 5=D 6=U)
	if (veto_id == 5) {
		double optical_coupling[13] = { 0., 0.94, 0.57, 0.35, 0.7, 0.094, 0.177, 0.52, 0.75, 0.52, 0.52, 0.38, 1.0 };
		for (int s = 0; s < 13; s++)
			optical_coupling[s] = optical_coupling[s] * 0.68;
		light_yield = 9200 / MeV;
		veff = 13 * cm / ns;
		sensor_effective_area = 0.9;
		sensor_qe = 0.25;
		sensor_gain = 1.;

		// Upper/lower
		if (channel == 1 || channel == 2) {
			// Get the paddle length: in veto paddles are along z
			length = aHit->GetDetector().dimensions[2];
			double s1 = aHit->GetDetector().dimensions[0];
			double s2 = aHit->GetDetector().dimensions[1];
			paddle_surface = 2 * s1 * 2 * s2;
			sensor_surface = pow(2.5 * cm, 2) * pi; // 2" pmt R-> (2*2.5/2 TBC)
			att_length = 350 * cm;
			light_guide_att = 1.0;
			// cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
		}
		if (channel == 5 || channel == 6) { // downstream upstream

			double s1 = aHit->GetDetector().dimensions[0];
			double s2 = aHit->GetDetector().dimensions[2];
			paddle_surface = 2 * s1 * 2 * s2; // surface perpendicular to the pmt position Surf=XxZ
			sensor_surface = pow((2.5 / 2) * cm, 2) * pi; //
			att_length = 400 * cm; // longer att lenght to take into account the perpendicular readout
			light_guide_att = 0.19; // no light guide
		}
		if (channel == 3 || channel == 4) {	//right left
		// Get the paddle length: in veto paddles are along y
			length = aHit->GetDetector().dimensions[1];
			double s1 = aHit->GetDetector().dimensions[0];
			double s2 = aHit->GetDetector().dimensions[2];
			sensor_surface = pow(2.5 * cm, 2) * pi; // 2" pmt R-> (2*2.5/2 TBC)
			paddle_surface = 2 * s1 * 2 * s2;
			att_length = 350 * cm;
			light_guide_att = 1.0;
			// cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
		}

		light_coll = sensor_surface / paddle_surface;
		if (sensor_surface > paddle_surface) light_coll = 1.;   // no more than the PMT size
		light_coll = optical_coupling[1] * light_coll * sensor_effective_area * light_guide_att; // coupling [1] identical for all

		// Get info about detector material to eveluate Birks effect
		double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
		//	cout << " Birks constant is: " << birks_constant << endl;
		//	cout << aHit->GetDetector().GetLogical()->GetMaterial()->GetName() << endl;

		double time_min[4] = { 0, 0, 0, 0 };

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();

		unsigned int nsteps = Edep.size();
		double Etot = 0;

		for (unsigned int s = 0; s < nsteps; s++)
			Etot = Etot + Edep[s];

		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double dLeft = -10000.;
				double dRight = -10000.;

				// Distances from left, right for upper/lower (along z)
				if (channel == 1 || channel == 2) {
					dLeft = length - Lpos[s].z();
					dRight = length + Lpos[s].z();
				}
				// Distances from top/bottom for side OV (along y)

				if (channel == 3 || channel == 4) {
					dLeft = length + Lpos[s].y();
					dRight = length - Lpos[s].y();
				}
				if (channel == 5 || channel == 6)
				// Distances from center for U/D OV (along y)

						{
					dLeft = Lpos[s].x();
					dRight = Lpos[s].y();
					double dCent = sqrt(dLeft * dLeft + dRight * dRight);
					dRight = dCent;

				}

				// cout << "\n Distances: " << endl;
				// cout << "\t dLeft, dRight " << dLeft <<  ", " << dRight << endl;

				// apply Birks effect
				// 			double stepl = 0.;

				//			if (s == 0){
				//				stepl = sqrt(pow((Lpos[s+1].x() - Lpos[s].x()),2) + pow((Lpos[s+1].y() - Lpos[s].y()),2) + pow((Lpos[s+1].z() - Lpos[s].z()),2));
				//			}
				//			else {
				//				stepl = sqrt(pow((Lpos[s].x() - Lpos[s-1].x()),2) + pow((Lpos[s].y() - Lpos[s-1].y()),2) + pow((Lpos[s].z() - Lpos[s-1].z()),2));
				//			}

				double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
				Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;
				// cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
				//	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;

				//if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with

				etotL = etotL + Edep_B / 2 * exp(-dLeft / att_length) * light_coll;
				etotR = etotR + Edep_B / 2 * exp(-dRight / att_length) * light_coll;

				//			  cout << "step: " << s << " etotL, etotR " << etotL << ", " << etotR  << endl;

				timeL = timeL + (times[s] + dLeft / veff) / nsteps;
				timeR = timeR + (times[s] + dRight / veff) / nsteps;

				if (etotL > 0.) {
					if (s == 0 || (time_min[0] > (times[s] + dLeft / veff))) time_min[0] = times[s] + dLeft / veff;
				}
				//      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
				if (etotR > 0.) {
					if (s == 0 || (time_min[1] > (times[s] + dRight / veff))) time_min[1] = times[s] + dRight / veff;
				}
			}
			//cout << " etotR " << etotR   <<  " ; etotL" <<  etotL <<endl;
			peL = G4Poisson(etotL * light_yield * sensor_qe);
			peR = G4Poisson(etotR * light_yield * sensor_qe);
			//cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
			//peL=(etotL*light_yield*sensor_qe);
			//peR=(etotR*light_yield*sensor_qe);
			//cout << " per " << peR   <<  " ; pel" <<  peL <<endl;

			double sigmaTL = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peL + 1.));
			double sigmaTR = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peR + 1.));
			//sigmaTL=0;
			//sigmaTR=0;
			tL = (time_min[0] + G4RandGauss::shoot(0.,sigmaTL))*1000.;                //time in ps
			tR = (time_min[1] + G4RandGauss::shoot(0.,sigmaTR))*1000.;                // time in ps
			// Digitization for ADC and QDC not used
			//TDC1=(int) (tL * tdc_conv);
			//TDC2=(int) (tR * tdc_conv);
			//if(TDC1<0) TDC1=0;
			//if(TDC2<0) TDC2=0;
			//ADC1=(int) (peL*sensor_gain*adc_conv + adc_ped);
			//ADC2=(int) (peR*sensor_gain*adc_conv + adc_ped);

			//	  cout << "ADC1: " << ADC1 << " " << peL << " " << sensor_gain << " " << adc_conv << endl;

			//cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
			//cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;

			//cout << " Light collection: " << light_coll << endl;

		}	// closes (Etot > 0) loop

		// ch 1,3               dRight
		// ch 2,4               dLeft
		// ch 7,8,9,10,11,12    dRight
		// ch 5,6               dRight
		// ignore the other side

		if (channel == 1 || channel == 2) {
			if (sector == 0) {
				ADC1 = peR;
				TDC1 = tR;
			} else if (sector == 1) {
				ADC1 = peL;
				TDC1 = tL;
			}
		}
		if (channel == 3 || channel == 4) {
			ADC1 = peR;
			TDC1 = tR;
		}
		if (channel == 5 || channel == 6) {
			ADC1 = peR;
			TDC1 = tR;
		}

		// cout << " ADC1: " << ADC1    <<  " ; TDC1: " <<  TDC1  << " ;  ADC2: "<< etot_g4*1000. << endl;
		if (verbosity > 4) {
			cout << log_msg << " veto: " << veto_id << ", channel: " << channel << ", sector: " << sector;
			cout << log_msg << " Etot=" << Etot / MeV << endl;
			cout << log_msg << " TDC1=" << TDC1 << " TDC2=" << TDC2 << " ADC1=" << ADC1 << " ADC2=" << ADC2 << endl;
		}

	} else if (veto_id == 4) {
		// proposal IV

		double veff = 13 * cm / ns;            // TO BE CHECKED
		// scintillator sizes
		double sx = aHit->GetDetector().dimensions[0];
		double sy = aHit->GetDetector().dimensions[1];
		double sz = aHit->GetDetector().dimensions[2];

//        double time_min[4] = {0,0,0,0};

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();
		unsigned int nsteps = Edep.size();
		double Etot = 0;
		double X_hit_ave = 0.;
		double Y_hit_ave = 0.;
		double Z_hit_ave = 0.;
		double T_hit_ave = 0.;
		double dLeft = -10000.;

		for (unsigned int s = 0; s < nsteps; s++)
			Etot = Etot + Edep[s];
		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;
				// average hit position XYZ
				X_hit_ave = X_hit_ave + Lpos[s].x();
				Y_hit_ave = Y_hit_ave + Lpos[s].y();
				Z_hit_ave = Z_hit_ave + Lpos[s].z();
				// average hit time
				T_hit_ave = T_hit_ave + times[s];

				//cout << "X " << Lpos[s].x() << " " << "Y " << Lpos[s].y() << " " << "Z " << Lpos[s].z() << " "<< "T " <<times[s] << " " << endl;

			}
			X_hit_ave = X_hit_ave / nsteps;
			Y_hit_ave = Y_hit_ave / nsteps;
			Z_hit_ave = Z_hit_ave / nsteps;
			T_hit_ave = T_hit_ave / nsteps;
			dLeft = sz - Z_hit_ave;
			timeL = dLeft / veff + T_hit_ave;

			double *pe_sipm;        // response for a mip (2..05 MeV energy released in 1cm thick)

			pe_sipm = IVresponseProposal(channel, X_hit_ave, Y_hit_ave, Z_hit_ave, sx, sy, sz);

			ADC1 = G4Poisson(pe_sipm[0] * etot_g4 / 2.05); // Scaling for more/less energy release)
			ADC2 = G4Poisson(pe_sipm[1] * etot_g4 / 2.05); // Scaling for more/less energy release)
			ADC3 = G4Poisson(pe_sipm[2] * etot_g4 / 2.05); // Scaling for more/less energy release)
			ADC4 = G4Poisson(pe_sipm[3] * etot_g4 / 2.05); // Scaling for more/less energy release)

			//adding a gaussian spread accoring to Luca's tabel
			ADC1 = (ADC1 + G4RandGauss::shoot(0.,13.));
			ADC2 = (ADC2 + G4RandGauss::shoot(0.,13.));
			ADC3 = (ADC3 + G4RandGauss::shoot(0.,13.));
			ADC4 = (ADC4 + G4RandGauss::shoot(0.,13.));
			if (ADC1 < 0) ADC1 = 0.;
			if (ADC2 < 0) ADC2 = 0.;
			if (ADC3 < 0) ADC3 = 0.;
			if (ADC4 < 0) ADC4 = 0.;

			double sigmaTL = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peL + 1.));
			sigmaTL = 0.;
			TDC1 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.; //time in ps
			TDC2 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.; //time in ps
			TDC3 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.; //time in ps
			TDC4 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.; //time in ps

			//cout <<  log_msg << " veto: " << veto_id   << ", channel: " << channel << ", sector: " << sector << endl;
			//cout << "X " << X_hit_ave << " " << "Y " << Y_hit_ave << " " << "Z " << Z_hit_ave << " "<< "T " <<T_hit_ave << " " << endl;
			//cout << "sipm1 " << pe_sipm[0] << " " << "sipm2 " << pe_sipm[1] << " " << "sipm3 " << pe_sipm[2] << " " << "sipm4 " << pe_sipm[3] << " " << endl;
			// cout << "dLeft " << dLeft << " " << "timeL " << timeL << " " << endl;

			//cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
			//cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;

			//cout << " Light collection: " << light_coll << endl;

		}
		// closes (Etot > 0) loop

	}

	// End Proposal

	// Outer veto PMT-readout (some id are missing since they are in a separte routine WLS readout)
	// 1 WLS readout (L/R)
	// 2 not existing anymore
	// 3 bottom
	else if (veto_id == 2 && channel != 1 && channel != 2 && channel != 5 && channel != 6) {
		//double optical_coupling[13]= {0., 0.94,0.57, 0.35, 0.7, 0.094, 0.177, 0.52, 0.75, 0.52, 0.52, 0.38, 1.0 };
		double optical_coupling[15] = {0., 0., 0., 0.35, 0.7, 0., 0., 0.94, 0.57, 0.52, 0.75, 0.52, 0.52, 0.38, 1.0};
		for (int s = 0; s < 15; s++)
		optical_coupling[s] = optical_coupling[s] * 0.68;
		light_yield = 9200 / MeV;
		veff = 13 * cm / ns;
		sensor_effective_area = 0.9;
		sensor_qe = 0.25;
		sensor_gain = 1.;

		// lower
		if (channel == 3 || channel == 4) {
			// Get the paddle length: in veto paddles are along z
			length = aHit->GetDetector().dimensions[2];
			double s1 = aHit->GetDetector().dimensions[0];
			double s2 = aHit->GetDetector().dimensions[1];
			paddle_surface = 2 * s1 * 2 * s2;
			sensor_surface = pow(2.5 * cm, 2) * pi;// 2" pmt R-> (2*2.5/2 TBC)
			att_length = 350 * cm;
			light_guide_att = 1.0;
			// cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
		}
		// 7-8-9-10 R, 11-12-13-14 L
		if (channel == 7 || channel == 8 || channel == 9 || channel == 10 || channel == 11 || channel == 12 || channel == 13 || channel == 14) {
			// Get the paddle length: in veto paddles are along y
			length = aHit->GetDetector().dimensions[1];
			double s1 = aHit->GetDetector().dimensions[0];
			double s2 = aHit->GetDetector().dimensions[2];
			sensor_surface = pow(2.5 * cm, 2) * pi;// 2" pmt R-> (2*2.5/2 TBC)
			paddle_surface = 2 * s1 * 2 * s2;
			att_length = 350 * cm;
			light_guide_att = 1.0;
			// cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
		}

		light_coll = sensor_surface / paddle_surface;
		if (sensor_surface > paddle_surface) light_coll = 1.;   // no more than the PMT size
		light_coll = optical_coupling[channel] * light_coll * sensor_effective_area * light_guide_att;// Including the coupling efficiency and the pc effective area
		//cout << " light collo " << light_coll     <<endl;

		// Get info about detector material to eveluate Birks effect
		double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
		//	cout << " Birks constant is: " << birks_constant << endl;
		//	cout << aHit->GetDetector().GetLogical()->GetMaterial()->GetName() << endl;

		double time_min[4] = {0, 0, 0, 0};

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();

		unsigned int nsteps = Edep.size();
		double Etot = 0;

		for (unsigned int s = 0; s < nsteps; s++)
		Etot = Etot + Edep[s];

		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double dLeft = -10000.;
				double dRight = -10000.;

				// Distances from left, right for upper/lower (along z)
				if (channel == 3 || channel == 4) {
					dLeft = length - Lpos[s].z();
					dRight = length + Lpos[s].z();
				}
				// Distances from left, right for other OV (along y)

				if (channel == 7 || channel == 8 || channel == 9 || channel == 10 || channel == 11 || channel == 12 || channel == 13 || channel == 14) {
					dLeft = length + Lpos[s].y();
					dRight = length - Lpos[s].y();
				}

				double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
				Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;

				etotL = etotL + Edep_B / 2 * exp(-dLeft / att_length) * light_coll;
				etotR = etotR + Edep_B / 2 * exp(-dRight / att_length) * light_coll;

				timeL = timeL + (times[s] + dLeft / veff) / nsteps;
				timeR = timeR + (times[s] + dRight / veff) / nsteps;

				if (etotL > 0.) {
					if (s == 0 || (time_min[0] > (times[s] + dLeft / veff))) time_min[0] = times[s] + dLeft / veff;
				}
				//      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
				if (etotR > 0.) {
					if (s == 0 || (time_min[1] > (times[s] + dRight / veff))) time_min[1] = times[s] + dRight / veff;
				}
			}
			//cout << " etotR " << etotR   <<  " ; etotL" <<  etotL <<endl;
			peL = G4Poisson(etotL * light_yield * sensor_qe);
			peR = G4Poisson(etotR * light_yield * sensor_qe);

			double sigmaTL = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peL + 1.));
			double sigmaTR = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peR + 1.));
			//sigmaTL=0;
			//sigmaTR=0;
			tL = (time_min[0] + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			tR = (time_min[1] + G4RandGauss::shoot(0.,sigmaTR))*1000.;// time in ps

		}
		// closes (Etot > 0) loop

		if (verbosity > 4) {
			cout << log_msg << " veto: " << veto_id << ", channel: " << channel;
			cout << log_msg << " Etot=" << Etot / MeV << endl;
			cout << log_msg << " TDC1=" << TDC1 << " TDC2=" << TDC2 << " ADC1=" << ADC1 << " ADC2=" << ADC2 << endl;
		}

		// ch 1,3               dRight
		// ch 2,4               dLeft
		// ch 7,8,9,10,11,12    dRight
		// ch 5,6               dRight
		// ignore the other side

		if (channel == 3) {
			ADC1 = peR;
			TDC1 = tR;
		}
		if (channel == 4) {
			ADC1 = peL;
			TDC1 = tL;
		}
		if (channel == 7 || channel == 8 || channel == 9 || channel == 10 || channel == 11 || channel == 12 || channel == 13 || channel == 14) {
			ADC1 = peR;
			TDC1 = tR;

		}
		// cout << " ADC1: " << ADC1    <<  " ; TDC1: " <<  TDC1  << " ;  ADC2: "<< etot_g4*1000. << endl;

	} // end of OV

	// Outer VETO OV WLS readout
	else if (veto_id == 2 && (channel == 1 || channel == 2 || channel == 5 || channel == 6)) {
		double veff = 13 * cm / ns; // TO BE CHECKED
		// scintillator sizes
		//        double sx=aHit->GetDetector().dimensions[0];
		//        double sy=aHit->GetDetector().dimensions[1];
		double sz = aHit->GetDetector().dimensions[2];

		//        double time_min[4] = {0,0,0,0};

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();
		unsigned int nsteps = Edep.size();
		double Etot = 0;
		double X_hit_ave = 0.;
		double Y_hit_ave = 0.;
		double Z_hit_ave = 0.;
		double T_hit_ave = 0.;
		double dLeft = -10000.;
		double dRight = -10000.;

		for (unsigned int s = 0; s < nsteps; s++)
		Etot = Etot + Edep[s];
		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;
				// average hit position XYZ
				X_hit_ave = X_hit_ave + Lpos[s].x();
				Y_hit_ave = Y_hit_ave + Lpos[s].y();
				Z_hit_ave = Z_hit_ave + Lpos[s].z();
				// average hit time
				T_hit_ave = T_hit_ave + times[s];

				//cout << "X " << Lpos[s].x() << " " << "Y " << Lpos[s].y() << " " << "Z " << Lpos[s].z() << " "<< "T " <<times[s] << " " << endl;

			}
			X_hit_ave = X_hit_ave / nsteps;
			Y_hit_ave = Y_hit_ave / nsteps;
			Z_hit_ave = Z_hit_ave / nsteps;
			T_hit_ave = T_hit_ave / nsteps;
			dLeft = sz - Z_hit_ave;
			dRight = sz + Z_hit_ave;
			timeL = dLeft / veff + T_hit_ave;
			timeR = dRight / veff + T_hit_ave;

			double *pe_wls;        // response for a mip (2..05 MeV energy released in 1cm thick)
			pe_wls = OVresponse(channel, X_hit_ave, Y_hit_ave, Z_hit_ave);

			ADC1 = G4Poisson(pe_wls[0] * etot_g4 / 5.0);// Scaling for more/less energy release 2.5cm = 5 MeV/MIPs )
			ADC2 = G4Poisson(pe_wls[1] * etot_g4 / 5.0);// Scaling for more/less energy release)
			//ADC3=G4Poisson(pe_wls[2]*etot_g4/5.0) ; // Scaling for more/less energy release)
			//ADC4=G4Poisson(pe_wls[3]*etot_g4/5.0) ; // Scaling for more/less energy release)
			//adding a gaussian spread  TBD
			//ADC1=(ADC1+G4RandGauss::shoot(0.,13.));
			//ADC2=(ADC2+G4RandGauss::shoot(0.,13.));
			//ADC3=(ADC3+G4RandGauss::shoot(0.,13.));
			//ADC4=(ADC4+G4RandGauss::shoot(0.,13.));
			if (ADC1 < 0) ADC1 = 0.;
			if (ADC2 < 0) ADC2 = 0.;
			if (ADC3 < 0) ADC3 = 0.;
			if (ADC4 < 0) ADC4 = 0.;

			double sigmaTL = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peL + 1.));
			sigmaTL = 0.;
			TDC1 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			TDC2 = (timeR + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			TDC3 = 0.;
			TDC4 = 0.;

			//   cout << "channel " << channel << endl;
			// cout << "X " << X_hit_ave << " " << "Y " << Y_hit_ave << " " << "Z " << Z_hit_ave << " "<< "T " <<T_hit_ave << " " << endl;
			//  cout << "sipm1 " << pe_sipm[0] << " " << "sipm2 " << pe_sipm[1] << " " << "sipm3 " << pe_sipm[2] << " " << "sipm4 " << pe_sipm[3] << " " << endl;
			// cout << "dLeft " << dLeft << " " << "timeL" << timeL << " " << endl;

			//cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
			//cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;

			//cout << " Light collection: " << light_coll << endl;

		}
		// closes (Etot > 0) loop
	}            // end of OV WLS readout

	// INNER VETO or CAL_PADs
	else if (veto_id == 1 || veto_id == 3) {            // || veto_id == 7 || veto_id == 8) {
		int chan = channel;
		if (veto_id == 3) {            //Cal_pad
			if (channel == 1) chan = 301;//top
			else if (channel == 2) chan = 302;//bottom
		}

		double veff = 13 * cm / ns;  // TO BE CHECKED
		// scintillator sizes
		double sx = aHit->GetDetector().dimensions[0];
		double sy = aHit->GetDetector().dimensions[1];
		double sz = aHit->GetDetector().dimensions[2];

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();
		unsigned int nsteps = Edep.size();
		double Etot = 0;
		double X_hit_ave = 0.;
		double Y_hit_ave = 0.;
		double Z_hit_ave = 0.;
		double T_hit_ave = 0.;
		double dLeft = -10000.;
		double dRight = -10000.;

		for (unsigned int s = 0; s < nsteps; s++) {
			Etot = Etot + Edep[s];
		}
		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;
				// average hit position XYZ
				X_hit_ave = X_hit_ave + Lpos[s].x();
				Y_hit_ave = Y_hit_ave + Lpos[s].y();
				Z_hit_ave = Z_hit_ave + Lpos[s].z();
				// average hit time
				T_hit_ave = T_hit_ave + times[s];
			}
			X_hit_ave = X_hit_ave / nsteps;
			Y_hit_ave = Y_hit_ave / nsteps;
			Z_hit_ave = Z_hit_ave / nsteps;
			T_hit_ave = T_hit_ave / nsteps;
			dLeft = sz - Z_hit_ave;
			dRight = sz + Z_hit_ave;
			timeL = dLeft / veff + T_hit_ave;
			timeR = dRight / veff + T_hit_ave;
			double *pe_sipm;		// response for a mip (2.05 MeV energy released in 1cm thick)
			pe_sipm = IVresponse(chan, X_hit_ave, Y_hit_ave, Z_hit_ave, sx, sy, sz);

			//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;

			ADC1 = G4Poisson(pe_sipm[0] * etot_g4 / 2.05);// Scaling for more/less ener54     gy release)
			ADC2 = G4Poisson(pe_sipm[1] * etot_g4 / 2.05);// Scaling for more/less energy release)
			ADC3 = G4Poisson(pe_sipm[2] * etot_g4 / 2.05);// Scaling for more/less energy release)
			ADC4 = G4Poisson(pe_sipm[3] * etot_g4 / 2.05);// Scaling for more/less energy release)
			//adding a gaussian spread accoring to Luca's tabel
			ADC1 = (ADC1 + G4RandGauss::shoot(0.,13.));
			ADC2 = (ADC2 + G4RandGauss::shoot(0.,13.));
			ADC3 = (ADC3 + G4RandGauss::shoot(0.,13.));
			ADC4 = (ADC4 + G4RandGauss::shoot(0.,13.));

			if (ADC1 < 0) ADC1 = 0.;
			if (ADC2 < 0) ADC2 = 0.;
			if (ADC3 < 0) ADC3 = 0.;
			if (ADC4 < 0) ADC4 = 0.;
			double sigmaTL = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peL + 1.));

			sigmaTL = 0.;
			TDC1 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			TDC2 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			TDC3 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			TDC4 = (timeL + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
		}
	}

	//BDX-HODO, BDX-MINI
	else if (veto_id == 6 || veto_id == 7 || veto_id == 8) {

		int chan = channel;
		if (veto_id == 6) chan = 600 + channel;			//bdx-hodo
		else if (veto_id == 7) chan = 700 + channel;//bdx-mini Outer Veto
		else if (veto_id == 8) chan = 800 + channel;//bdx-mini Inner Veto

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();
		unsigned int nsteps = Edep.size();

		double Etot = 0;
		double X_hit_ave = 0.;
		double Y_hit_ave = 0.;
		double Z_hit_ave = 0.;
		double T_hit_ave = 0.;
		for (unsigned int s = 0; s < nsteps; s++) {
			Etot = Etot + Edep[s];
		}
		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;
				// average hit time
				T_hit_ave = T_hit_ave + times[s];
			}
			T_hit_ave = T_hit_ave / nsteps;

			// Different procedure for BDX-Hodo: using directly the extracted pe from ps_sipm
			if (chan >= 600 & chan < 700) {
				//double MPV[14]={0.,   205.3, 155.6, 179.4, 177.4, 198.8, 184.9, 208.5, 207.1, 226.9, 213.1, 227.4, 164.7, 208.1 };
				double MPV[14] = {0., 205.3, 155.6, 179.4, 177.4, 198.8, 184.9, 208.5, 207.1, 226.9, 213.1, 227.4, 164.7, 208.1};
				double S_Land[14] = {0., 9.7, 10.9, 11.5, 5.0, 10.2, 9.3, 4.5, 7.7, 4.3, 4.6, 5.2, 7.1, 3.8};
				double S_Gaus[14] = {0., 19.3, 16.1, 20.2, 37.0, 21.1, 17.8, 27.5, 19.8, 21.5, 24.6, 15.2, 23.0, 24.5};
				ADC1 = G4Poisson(1.2 * MPV[channel] * etot_g4 / 2.05); //this is not correct for Ch=12 and 13 being thickness=2cm (and not 1cm)
				//ADC1=(ADC1+G4RandGauss::shoot(0.,(S_Land[channel]+S_Gaus[channel])));
				ADC1 = (ADC1 + G4RandGauss::shoot(0.,10.));
				ADC2 = 0.;
				ADC3 = 0.;
				ADC4 = 0.;
				//cout <<  " ++ HIT BEGIN ++++++" << endl ;
				//cout <<  " chan: " << channel << endl ;
				//cout <<  " ADC1: " << ADC1 << endl ;
				//cout <<  " ++ HIT END ++++++" << channel << endl ;
			}

			//For BDX-MINI cylinder/octagon vetos (chan==701 || chan ==801) sum up lights in the position of the 8 sipms step by step
			else if (chan >= 700 & chan < 900) { //BDX-MINI vetos
				if (chan == 701 || chan == 801) { //BDX-MINI Cylinder/Octagon
					//finding time clusters
					auto ghits=sortBDXMiniVetoHits(aHit);

					//  cout <<" "<< endl;
					double QSipmBdxMini[8];
					double TSipmBdxMini[8];

					double att_l_ang;
					double att_z;
					double TrGrv;
					//	double *MeV2pe;
					double MeV2pe[8];

					double Qdep;
					unsigned int NGrv;
					double PhiLoc;
					double SigmaTSipm = 0.;// Time spread on sipm on octagon and cylinder
					double QTThre = 10.;//Timing is considered only if the j-mo sipm receive a charge larger than QTThre pe
					for (unsigned int s = 0; s < 8; s++) {
						QSipmBdxMini[s] = 0.;
						TSipmBdxMini[s] = 1000.;                       // Time initializatiom
					}

					//looping on found clusters
					for (auto ghit:ghits) {

						double phiCluster=(atan2(ghit.avgLPos.x(),ghit.avgLPos.y()) / acos(-1.) * 180. + 180);
						double zCluster=ghit.avgLPos.z();
						double DeltaZ=0;
						//  cout << "CLUSTER  N=" <<s<< endl;
						for (unsigned int j = 0; j < 8; j++) {                       //looping on 8 sipm
							if (chan == 701) { //OV                      // outer veto sipm staggered by 22.5deg
								PhiLoc = 360 - abs((j) * 45. - phiCluster - 22.5);
								if (PhiLoc < 0) PhiLoc = 360 + PhiLoc;
								//	MeV2pe = MeV2peOV;
								MeV2pe[j] = MeV2peOV[j]*tunOV[j];
								TrGrv = TrGrvOL;
								att_l_ang = att_l_ang_OV;
								att_z = att_z_OV;
								DeltaZ=zCluster+aHit->GetDetector().dimensions[2];
							}
							else if (chan == 801) { //IV

								auto PhiSipm = -22.5 + j*45;//A.C added here 22.5 due to the rotation in geometry
								PhiLoc = fabs(phiCluster - PhiSipm);
								if (PhiLoc>360) PhiLoc=PhiLoc-360;
								if (PhiLoc<0) PhiLoc=PhiLoc+360;
								MeV2pe[j] = MeV2peIV[j]*tunIV[j];
								TrGrv = TrGrvIL;
								att_l_ang = att_l_ang_IV;
								att_z = att_z_IV;
								DeltaZ=zCluster+aHit->GetDetector().dimensions[9];
							}
							NGrv = int(PhiLoc / 45.); // counting grooves
							//A.C. adding z-dependence
							//PhiLoc: angular difference between SiPM and cluster in degrees
							//DeltaZ: longitudal distance between cluster and SiPms bottom plane in cm
							double attZ_L=1;
							double attZ_R=1;
							if (PhiLoc<180) {
								attZ_L=1-(1-exp(-DeltaZ/att_z))*cos(PhiLoc*3.14156/(2*180.));
							}
							if (PhiLoc>180) {
								attZ_R=1-(1-exp(-DeltaZ/att_z))*cos((360.-PhiLoc)*3.14156/(2*180.));
							}
							//	Qdep = MeV2pe[j] * ghit.Esum * (pow(TrGrv, NGrv) * exp(-PhiLoc / att_l_ang) * attZ_L + pow(TrGrv, (7 - NGrv)) * exp(-abs(360. - PhiLoc) / att_l_ang) * attZ_R);
							//TEST ATT//
							Qdep = MeV2pe[j] * ghit.Esum * (pow(TrGrv, NGrv) * exp(-PhiLoc / att_l_ang) * attZ_L + pow(TrGrv, (7 - NGrv)) * exp(-abs(360. - PhiLoc) / att_l_ang) * attZ_R);


							QSipmBdxMini[j] = QSipmBdxMini[j] + Qdep;


							// Timing
							double DeltaPhi = abs(PhiLoc);
							if (Qdep > QTThre) { // Timing is considered only if the j-mo sipm receive a charge larger than QTThre pe
								double TOld;
								double TNew;
								if (DeltaPhi > (360 - PhiLoc)) DeltaPhi = abs(360 - PhiLoc);
								TOld = TSipmBdxMini[j];
								TNew = ghit.avgT + abs(DeltaPhi / LightSpeedAng);
								if (TNew < TOld && TNew < TSipmBdxMini[j]) TSipmBdxMini[j] = TNew;
							}
						}
					}

					for (unsigned int s = 0; s < 8; s++) {
						QSipmBdxMini[s] = G4Poisson(QSipmBdxMini[s]);
					}

					for (unsigned int s = 0; s < 8; s++) {
						TSipmBdxMini[s] = (TSipmBdxMini[s] + G4RandGauss::shoot(0.,SigmaTSipm))*1000.;                      //time in ps
					}

					ADC1 = QSipmBdxMini[0];
					ADC2 = QSipmBdxMini[1];
					ADC3 = QSipmBdxMini[2];
					ADC4 = QSipmBdxMini[3];
					ADC5 = QSipmBdxMini[4];
					ADC6 = QSipmBdxMini[5];
					ADC7 = QSipmBdxMini[6];
					ADC8 = QSipmBdxMini[7];

					TDC1 = TSipmBdxMini[0];
					TDC2 = TSipmBdxMini[1];
					TDC3 = TSipmBdxMini[2];
					TDC4 = TSipmBdxMini[3];
					TDC5 = TSipmBdxMini[4];
					TDC6 = TSipmBdxMini[5];
					TDC7 = TSipmBdxMini[6];
					TDC8 = TSipmBdxMini[7];
				}
				/*Checked by A.C.: 709->OV TOP, 809->IV TOP, 810->IV BOTTOM, 710->OV BOTTOM*/
				else if (chan == 709 || chan == 710 || chan == 809 || chan == 810) {
					// BDX-MINI vetos (OuterTop (RECON: channel-10), OuterBottom (RECON: channel-9), InnerTop(RECON: channel-10), InnerBottom(RECON: channel-9)) leads
					double LY[4] = {83., 76., 112.5, 154.5};
					double QSipmLeadsBdxMini = 0.;
					double TSipmLeadsBdxMini = 0.;
					unsigned int j = 0;
					double SigmaTSipmLeads = 0.;		// Time spread on sipm leads
					if (chan == 709) j = 0;//OV-TOP
					if (chan == 710) j = 1;//OV-BOTTOM
					if (chan == 809) j = 2;//IV-TOP
					if (chan == 810) j = 3;//IV-BOTTOM

					QSipmLeadsBdxMini = G4Poisson(LY[j] * etot_g4);
					double sigmaTL = 0.;
					TSipmLeadsBdxMini = (T_hit_ave + G4RandGauss::shoot(0.,SigmaTSipmLeads))*1000.;//time in ps

					//         cout << "ADC leads " << chan<< " "<< QSipmLeadsBdxMini << " " << "TDC leads " << TSipmLeadsBdxMini <<endl;
					ADC1 = QSipmLeadsBdxMini;
					TDC1 = TSipmLeadsBdxMini;
				}
			}	// End BDX-MINI veto response
		}	// closes (Etot > 0) loop
	}
//starting paddles
	else if (veto_id == 4) {
		double optical_coupling[3] = {0., 1., 0.37};
		for (int s = 0; s < 3; s++)
		optical_coupling[s] = optical_coupling[s] * 0.34;

		light_yield = 9200 / MeV;
		veff = 13 * cm / ns;
		sensor_surface = pow(1.27 * cm, 2) * pi; // 1" pmt R-> (2.5/2 TBC)
		sensor_effective_area = 0.9;

		sensor_qe = 0.25;
		sensor_gain = 1.;

		// Get the paddle length: in veto paddles are along z
		length = aHit->GetDetector().dimensions[2];
		double s1 = aHit->GetDetector().dimensions[0];
		double s2 = aHit->GetDetector().dimensions[1];
		paddle_surface = 2 * s1 * 2 * s2;
		att_length = 350 * cm;
		light_guide_att = 1.;
		//cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;

		light_coll = sensor_surface / paddle_surface;
		if (sensor_surface > paddle_surface) light_coll = 1.;// no more than the PMT size
		light_coll = optical_coupling[channel] * light_coll * sensor_effective_area * light_guide_att;// Including the coupling efficiency and the pc effective area
		// cout << " channel,veto: " << channel << " " << veto_id    <<  " optical-couping: " <<  optical_coupling[channel]  <<  " light coll: " <<  light_coll  <<endl;
		//cout << " light collo " << light_coll     <<endl;

		// Get info about detector material to eveluate Birks effect
		double birks_constant = aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
		//	cout << " Birks constant is: " << birks_constant << endl;
		//	cout << aHit->GetDetector().GetLogical()->GetMaterial()->GetName() << endl;

		double time_min[4] = {0, 0, 0, 0};

		vector<G4ThreeVector> Lpos = aHit->GetLPos();
		vector<G4double> Edep = aHit->GetEdep();
		vector<G4double> Dx = aHit->GetDx();
		// Charge for each step
		vector<int> charge = aHit->GetCharges();
		vector<G4double> times = aHit->GetTime();

		unsigned int nsteps = Edep.size();
		double Etot = 0;

		for (unsigned int s = 0; s < nsteps; s++)
		Etot = Etot + Edep[s];
		//for(unsigned int s=0; s<nsteps; s++) cout << "Energy = " << Edep[s]*1000.*1000. << endl;

		if (Etot > 0) {
			for (unsigned int s = 0; s < nsteps; s++) {
				double dLeft = -10000.;
				double dRight = -10000.;

				// Distances from left, right for upper/lower (along z)
				dLeft = length - Lpos[s].z();
				dRight = length + Lpos[s].z();

				//cout << "\n Distances: " << endl;
				// cout << "\t dLeft, dRight " << dLeft <<  ", " << dRight << endl;

				// apply Birks effect
				// 			double stepl = 0.;

				//			if (s == 0){
				//				stepl = sqrt(pow((Lpos[s+1].x() - Lpos[s].x()),2) + pow((Lpos[s+1].y() - Lpos[s].y()),2) + pow((Lpos[s+1].z() - Lpos[s].z()),2));
				//			}
				//			else {
				//				stepl = sqrt(pow((Lpos[s].x() - Lpos[s-1].x()),2) + pow((Lpos[s].y() - Lpos[s-1].y()),2) + pow((Lpos[s].z() - Lpos[s-1].z()),2));
				//			}

				double Edep_B = BirksAttenuation(Edep[s], Dx[s], charge[s], birks_constant);
				Edep_B = Edep[s];
				etot_g4 = etot_g4 + Edep_B;
				// cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
				//	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;

				//if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with

				etotL = etotL + Edep_B / 2 * exp(-dLeft / att_length) * light_coll;
				etotR = etotR + Edep_B / 2 * exp(-dRight / att_length) * light_coll;

				//			  cout << "step: " << s << " etotL, etotR " << etotL << ", " << etotR  << endl;

				timeL = timeL + (times[s] + dLeft / veff) / nsteps;
				timeR = timeR + (times[s] + dRight / veff) / nsteps;

				if (etotL > 0.) {
					if (s == 0 || (time_min[0] > (times[s] + dLeft / veff))) time_min[0] = times[s] + dLeft / veff;
				}
				//      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
				if (etotR > 0.) {
					if (s == 0 || (time_min[1] > (times[s] + dRight / veff))) time_min[1] = times[s] + dRight / veff;
				}
			}
			//cout << " etotR " << etotR   <<  " ; etotL" <<  etotL <<endl;
			peL = G4Poisson(etotL * light_yield * sensor_qe);
			peR = G4Poisson(etotR * light_yield * sensor_qe);
			//cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
			//peL=(etotL*light_yield*sensor_qe);
			//peR=(etotR*light_yield*sensor_qe);
			//cout << " per " << peR   <<  " ; pel" <<  peL <<endl;

			double sigmaTL = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peL + 1.));
			double sigmaTR = sqrt(pow(0.2 * nanosecond, 2.) + pow(1. * nanosecond, 2.) / (peR + 1.));
//            sigmaTL=0;
//            sigmaTR=0;
			tL = (time_min[0] + G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
			tR = (time_min[1] + G4RandGauss::shoot(0.,sigmaTR))*1000.;// time in ps
			//  cout << " tL " << tL   <<  " ; timeL" <<  timeL <<endl;
			// Digitization for ADC and QDC not used
			//TDC1=(int) (tL * tdc_conv);
			//TDC2=(int) (tR * tdc_conv);
			//if(TDC1<0) TDC1=0;
			//if(TDC2<0) TDC2=0;
			//ADC1=(int) (peL*sensor_gain*adc_conv + adc_ped);
			//ADC2=(int) (peR*sensor_gain*adc_conv + adc_ped);

			//	  cout << "ADC1: " << ADC1 << " " << peL << " " << sensor_gain << " " << adc_conv << endl;

			//cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
			//cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;

			//cout << " Light collection: " << light_coll << endl;

		}
		// closes (Etot > 0) loop

		if (verbosity > 4) {
			cout << log_msg << " veto: " << veto_id << ", channel: " << channel;
			cout << log_msg << " Etot=" << Etot / MeV << endl;
			cout << log_msg << " TDC1=" << TDC1 << " TDC2=" << TDC2 << " ADC1=" << ADC1 << " ADC2=" << ADC2 << endl;
		}

		// ch 1,3               dRight
		// ch 2,4               dLeft
		// ch 7,8,9,10,11,12    dRight
		// ch 5,6               dRight
		// ignore the other side

		ADC1 = peR;
		TDC1 = tR;
		// cout << " ADC1: " << ADC1    <<  " ; TDC1: " <<  TDC1  << " ;  ADC2: "<< etot_g4*1000. << endl;

	}        // end of paddles

	dgtz["hitn"] = hitn;
	dgtz["sector"] = sector;
	dgtz["veto"] = veto_id;
	dgtz["channel"] = channel;
	dgtz["adc1"] = ADC1;        // output in pe
	dgtz["adc2"] = ADC2;        //deposited energy in keV
	dgtz["adc3"] = ADC3;        // ignore
	dgtz["adc4"] = ADC4;        // ignore
	dgtz["adc5"] = ADC5;        // ignore
	dgtz["adc6"] = ADC6;        // ignore
	dgtz["adc7"] = ADC7;        // ignore
	dgtz["adc8"] = ADC8;        // ignore

	dgtz["tdc1"] = TDC1;        // output in ps
	dgtz["tdc2"] = TDC2;        // ignore
	dgtz["tdc3"] = TDC3;        // ignore
	dgtz["tdc4"] = TDC4;        // ignore
	dgtz["tdc5"] = TDC5;        // ignore
	dgtz["tdc6"] = TDC6;        // ignore
	dgtz["tdc7"] = TDC7;        // ignore
	dgtz["tdc8"] = TDC8;        // ignore

	return dgtz;
}

vector<identifier> veto_HitProcess::processID(vector<identifier> id, G4Step *step, detector Detector) {
	id[id.size() - 1].id_sharing = 1;
	return id;
}

double* veto_HitProcess::IVresponse(int channel, double xx, double yy, double zz, double sx, double sy, double sz) {

// Response of the different IV plastic paddles
// ch
//
//
	static double response[4];

	for (unsigned int s = 0; s < 4; s++)
		response[s] = 0.;

	if (channel == 1000)        //top # Run 1 - fall 2015 - fall 2016 facking channel=1000
			{
		double x = xx / 10.;
		double y = (sz + zz) / 10.;
		double normfactor[4] = { 1.25, 1.65, 1.26, 2.21 };

		double parm[4][8] = { { 1.99627e+01, 1.64910e-01, -5.83528e-01, -7.34483e-03, -1.25062e-03, 4.43805e-03, 5.63766e-05, 1.40682e-05 }, { 1.86162e+01, 4.36475e-02, -6.78752e-02, -5.47887e-03, -1.60512e-04, -2.33958e-02, 5.55285e-05, -5.94424e-05 }, { 1.85966e+01,
				1.96301e-01, 1.34868e-01, -7.66131e-04, -1.61720e-03, -1.91598e-02, -1.76198e-06, -4.72970e-05 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}
		//     cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 2000)        //bottom # Run 1 - fall 2015 - fall 2016 facking channel=2000
			{        // Assuming an overall size of 42.8 cm with 4 bars of
		double x = -(xx - 428 / 2) / 10;
		double y = (sz + zz) / 10.;
		double normfactor[4] = { 1.55, 3.9, 2.92, 2.75 };

		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		if (x < 10) response[0] = (-0.000303034) * y * y + (0.00658939) * y + 32.4847; //D1
		if (x > 10 && x < 20) response[1] = (0.00301674) * y * y + (-0.446544) * y + 27.6374; //D4
		if (x > 20 && x < 32.8) response[2] = (-0.000275694) * y * y + (0.00124251) * y + 18.8999; //D3
		if (x > 32.8 && x < 42.8) response[3] = (-0.00139525) * y * y + (0.104993) * y + 18.1047; //D2
		for (unsigned int s = 0; s < 4; s++)
			response[s] = response[s] * normfactor[s];

		// cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	}

	if (channel == 1)      //top Dec16: waiting for the new pmap I'm using the old top (now positioned as bottom)
			{
		double x = xx / 10.;
		double y = (sz - zz) / 10.;
		double normfactor[4] = { 1.25, 1.65, 1.26, 2.21 };

		double parm[4][8] = { { 1.99627e+01, 1.64910e-01, -5.83528e-01, -7.34483e-03, -1.25062e-03, 4.43805e-03, 5.63766e-05, 1.40682e-05 }, { 1.86162e+01, 4.36475e-02, -6.78752e-02, -5.47887e-03, -1.60512e-04, -2.33958e-02, 5.55285e-05, -5.94424e-05 }, { 1.85966e+01,
				1.96301e-01, 1.34868e-01, -7.66131e-04, -1.61720e-03, -1.91598e-02, -1.76198e-06, -4.72970e-05 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}

		//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;
		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 2)        //bottom # Run 2 - fall 2016 - current | This is the Run 1 top now positioned on the bottom (assuming we moved rigidly down the paddle)
			{
		double x = xx / 10.;
		double y = (sz - zz) / 10.;
		double normfactor[4] = { 1.25, 1.65, 1.26, 2.21 };

		double parm[4][8] = { { 1.99627e+01, 1.64910e-01, -5.83528e-01, -7.34483e-03, -1.25062e-03, 4.43805e-03, 5.63766e-05, 1.40682e-05 }, { 1.86162e+01, 4.36475e-02, -6.78752e-02, -5.47887e-03, -1.60512e-04, -2.33958e-02, 5.55285e-05, -5.94424e-05 }, { 1.85966e+01,
				1.96301e-01, 1.34868e-01, -7.66131e-04, -1.61720e-03, -1.91598e-02, -1.76198e-06, -4.72970e-05 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}

		//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;
		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 3)        // Side Upstream
			{        // Assuming an overall size of 42.8 cm with 4 bars of
		double x = -xx / 10;
		double y = (yy + sy) / 10.;

		double normfactor[1] = { 1.13 };

		double parm[4] = { -0.04, -0.05, 1.4, 85. };

		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		response[0] = parm[0] * x * x + parm[1] * y * y + parm[2] * y + parm[3];
		for (unsigned int s = 0; s < 1; s++)
			response[s] = response[s] * normfactor[s];

		//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;
		//cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 4)        // Side Downstream
			{        // Assuming an overall size of 42.8 cm with 4 bars of
		double x = -xx / 10;
		double y = (yy + sy) / 10.;

		double parm[4] = { -0.04, -0.05, 1.4, 75. };
		double normfactor[1] = { 1.03 };

		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		response[0] = parm[0] * x * x + parm[1] * y * y + parm[2] * y + parm[3];
		for (unsigned int s = 0; s < 1; s++)
			response[s] = response[s] * normfactor[s];

		//cout << "Downstream" << endl;
		//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;
		//cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	}

	else if (channel == 5)        //Right
			{
		double x = -yy / 10.;
		double y = (sz - zz) / 10.;

		double parm[4][8] = { { 2.34524e+01, 4.28317e-02, -5.91894e-01, -5.13309e-03, -2.47905e-04, -3.44887e-03, 4.25481e-05, -1.03817e-05 }, { 1.68313e+01, 5.36853e-02, -2.14037e-01, -4.80535e-03, -4.65364e-04, -1.66572e-02, 4.89028e-05, -3.33380e-05 }, { 2.50310e+01,
				-3.10007e-02, 3.57657e-01, -1.39833e-02, 2.99406e-04, -3.23669e-02, 1.27237e-04, -1.13100e-05 }, { 1.74834e+01, 1.83925e-01, 5.36737e-01, 7.09769e-04, -1.64490e-03, 7.48199e-03, 3.43011e-08, 2.11894e-05 }

		};
		double normfactor[4] = { 0.97, 1.25, 1.06, 1.09 };

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}

		//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;
		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 6)        //Left
			{
		double x = -yy / 10.;
		double y = (sz - zz) / 10.;

		double parm[4][8] = { { 8.12418e+00, 6.61315e-02, -2.99641e-01, -9.10408e-04, -6.79474e-04, 2.00648e-03, 1.24963e-05, -1.73809e-05 }, { 1.19501e+01, 4.76291e-02, -1.77047e-01, 9.27111e-05, -4.63061e-04, -1.40014e-02, 4.39766e-06, -2.93896e-05 }, { 1.68607e+01,
				-4.15476e-02, 2.54857e-01, -6.87363e-03, 3.26876e-04, -2.65178e-02, 5.62748e-05, -3.56067e-06 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};
		double normfactor[4] = { 1.27, 1.14, 1.0, 0.83 };

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}
		//cout << "Left" << endl;
		//cout << "sx " << sx << " " << "sy " << sy << " " << "sz " << sz << endl;
		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 301 || channel == 302)        // Cal_pad top/bottom TB refined with measured response. So far using old parametrization of D1
			{        // Assuming an overall size of 10cm x 40cm
					 //double x=-(xx-sx)/10;
		double y = (sz - zz) / 10.;
		double normfactor[4] = { 1.55, 3.9, 2.92, 2.75 };

		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		response[0] = (-0.000303034) * y * y + (0.00658939) * y + 32.4847; //D1
		//if (x>10 && x <20 ) response[1]=   (0.00301674)*y*y + (-0.446544)*y + 27.6374; //D4
		//if (x>20 && x <32.8 ) response[2]= (-0.000275694)*y*y + (0.00124251)*y + 18.8999; //D3
		//if (x>32.8 && x <42.8 ) response[3]= (-0.00139525)*y*y + (0.104993)*y + 18.1047; //D2
		for (unsigned int s = 0; s < 4; s++)
			response[s] = response[s] * normfactor[s];

		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel >= 600 & channel < 700)        // BDX-hodo NOT USED
			{        // Assuming an overall size of 10cm x 40cm
		double x = -(xx - sx) / 10;
		double y = (sz - zz) / 10.;
		// BDX-Hodo scint parameters
		//double MPV[14]={0.,  205.3, 155.6, 179.4, 177.4, 198.8, 184.9, 208.5, 207.1, 226.9, 213.1, 227.4, 164.7, 208.1 };
		//double S_Land[14]={0., 9.7,  10.9,  11.5,   5.0,  10.2,   9.3,   4.5,   7.7,   4.3,   4.6,   5.2,   7.1,  3.8};
		//double S_Gaus[14]={0.  19.3, 16.1,  20.2,  37.0, 2 1.1,  17.8,  27.5,  19.8,  21.5,  24.6,  15.2,  23.0, 24.5};
		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		double normfactor[4] = { 1.55, 3.9, 2.92, 2.75 };
		response[0] = 0.;

		for (unsigned int s = 0; s < 4; s++)
			response[s] = response[s] * normfactor[s];
		//cout <<  " ++ HIT BEGIN ++++++" << endl ;
		//cout <<  " chan: " << channel << endl ;
		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
		//cout <<  " ++ HIT END ++++++" << channel << endl ;

	}
//  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
// cout <<  " res[0]: " << response[0] << " channel " << channel<<endl;
	return response;
}

double* veto_HitProcess::IVresponseProposal(int channel, double xx, double yy, double zz, double sx, double sy, double sz) {
// Response of the different IV plastic paddles
// ch
//
//
	static double response[4];

	for (unsigned int s = 0; s < 4; s++)
		response[s] = 0.;

	if (channel == 1)			//top
			{
		double x = xx / 10.;
		double y = (sz - zz) / 10.;
		double normfactor[4] = { 1.23, 1.45, 1.25, 1.91 };

		double parm[4][8] = { { 1.99627e+01, 1.64910e-01, -5.83528e-01, -7.34483e-03, -1.25062e-03, 4.43805e-03, 5.63766e-05, 1.40682e-05 }, { 1.86162e+01, 4.36475e-02, -6.78752e-02, -5.47887e-03, -1.60512e-04, -2.33958e-02, 5.55285e-05, -5.94424e-05 }, { 1.85966e+01,
				1.96301e-01, 1.34868e-01, -7.66131e-04, -1.61720e-03, -1.91598e-02, -1.76198e-06, -4.72970e-05 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}
		//     cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 2)        //bottom copying the top with X,Y swapped
			{        // Assuming an overall size of 42.8 cm with 4 bars of
		double x = -(xx - sx) / 10;
		double y = (sz - zz) / 10.;

		double normfactor[4] = { 1.23, 1.45, 1.25, 1.91 };

		double parm[4][8] = { { 1.99627e+01, 1.64910e-01, -5.83528e-01, -7.34483e-03, -1.25062e-03, 4.43805e-03, 5.63766e-05, 1.40682e-05 }, { 1.86162e+01, 4.36475e-02, -6.78752e-02, -5.47887e-03, -1.60512e-04, -2.33958e-02, 5.55285e-05, -5.94424e-05 }, { 1.85966e+01,
				1.96301e-01, 1.34868e-01, -7.66131e-04, -1.61720e-03, -1.91598e-02, -1.76198e-06, -4.72970e-05 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}

		// cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 3)        // Side Upstream
			{        // Assuming an overall size of 42.8 cm with 4 bars of
		double x = -xx / 10;
		double y = (yy + sy) / 10.;

		double parm[4] = { -0.04, -0.05, 1.4, 85. };
		double normfactor[1] = { 1.13 };

		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		response[0] = parm[0] * x * x + parm[1] * y * y + parm[2] * y + parm[3];
		for (unsigned int s = 0; s < 1; s++)
			response[s] = response[s] * normfactor[s];

		// cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
		// cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 4)        // Side Downstream
			{        // Assuming an overall size of 42.8 cm with 4 bars of
		double x = xx / 10;
		double y = (yy + sy) / 10.;

		double parm[4] = { -0.04, -0.05, 1.4, 75. };
		double normfactor[1] = { 1.03 };

		for (unsigned int s = 0; s < 4; s++)
			response[s] = 0.;
		response[0] = parm[0] * x * x + parm[1] * y * y + parm[2] * y + parm[3];
		for (unsigned int s = 0; s < 1; s++)
			response[s] = response[s] * normfactor[s];

		//cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	}

	else if (channel == 5)        //Right
			{
		double x = -yy / 10.;
		double y = (sz - zz) / 10.;

		double parm[4][8] = { { 2.34524e+01, 4.28317e-02, -5.91894e-01, -5.13309e-03, -2.47905e-04, -3.44887e-03, 4.25481e-05, -1.03817e-05 }, { 1.68313e+01, 5.36853e-02, -2.14037e-01, -4.80535e-03, -4.65364e-04, -1.66572e-02, 4.89028e-05, -3.33380e-05 }, { 2.50310e+01,
				-3.10007e-02, 3.57657e-01, -1.39833e-02, 2.99406e-04, -3.23669e-02, 1.27237e-04, -1.13100e-05 }, { 1.74834e+01, 1.83925e-01, 5.36737e-01, 7.09769e-04, -1.64490e-03, 7.48199e-03, 3.43011e-08, 2.11894e-05 }

		};
		double normfactor[4] = { 1.08, 1.02, 0.96, 0.97 };

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}
		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

	} else if (channel == 6)        //Left
			{
		double x = -yy / 10.;
		double y = (sz - zz) / 10.;

		double parm[4][8] = { { 8.12418e+00, 6.61315e-02, -2.99641e-01, -9.10408e-04, -6.79474e-04, 2.00648e-03, 1.24963e-05, -1.73809e-05 }, { 1.19501e+01, 4.76291e-02, -1.77047e-01, 9.27111e-05, -4.63061e-04, -1.40014e-02, 4.39766e-06, -2.93896e-05 }, { 1.68607e+01,
				-4.15476e-02, 2.54857e-01, -6.87363e-03, 3.26876e-04, -2.65178e-02, 5.62748e-05, -3.56067e-06 }, { 9.73394e+00, 1.56111e-01, 3.27558e-01, 2.45041e-03, -1.31615e-03, 5.82688e-03, -1.48528e-05, 2.35177e-05 }

		};
		double normfactor[4] = { 1.09, 1.14, 0.68, 0.63 };

		for (unsigned int s = 0; s < 4; s++) {
			response[s] = parm[s][7] * x * x * y + parm[s][6] * x * y * y + parm[s][5] * x * x + parm[s][4] * y * y + parm[s][3] * x * y + parm[s][2] * x + parm[s][1] * y + parm[s][0];
			response[s] = response[s] * normfactor[s];
		}

	}

//  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
// cout <<  " res[0]: " << response[0] << " channel " << channel<<endl;
	return response;
}

double* veto_HitProcess::OVresponse(int channel, double xx, double yy, double zz) {

// Response of the different OV plastic paddles read by WLS + PMT
// ch
//
//
	static double response[4];

	for (unsigned int s = 0; s < 4; s++)
		response[s] = 0.;

	if (channel == 1)        // TOp read by Left and Right sides
			{
		double pd_len = 1810 / 10.; //in cm
		//double x=xx/10.;
		double y = pd_len / 2 + zz / 10;			//in cm
		double normfactor[2] = { 1.0, 1.0 };

		double parm[2][2] = { { 178., 177.3 }, { 200., 259 } };

		response[0] = normfactor[0] * parm[0][0] * exp(-y / parm[0][1]); // Right/Upstream
		response[1] = normfactor[1] * parm[1][0] * exp(-(pd_len - y) / parm[1][1]); // Left/Downstream

		//    cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
	}
	if (channel == 5)        // Upstream p read by Left and Right sides
			{
		double pd_len = 668 / 10.; //in cm
		//double x=xx/10.;
		double y = pd_len / 2 + yy / 10;			//in cm
		double normfactor[2] = { 1.0, 1.0 };

		double parm[2][2] = { { 245., 480 }, { 245., 480 } };

		response[0] = normfactor[0] * parm[0][0] * exp(-y / parm[0][1]); // Right/bottom
		response[0] = 0.; // Only TOP readout
		response[1] = normfactor[1] * parm[1][0] * exp(-(pd_len - y) / parm[1][1]); // Left/top
		response[0] = response[1]; // Forcing ADC1 to have the only readout channel
		response[1] = 0.;

		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		//cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
	}
	if (channel == 6)        // Upstream p read by Left and Right sides
			{
		double pd_len = 668 / 10.; //in cm
		//double x=xx/10.;
		double y = pd_len / 2 + yy / 10;			//in cm
		double normfactor[2] = { 1.0, 1.0 };

		double parm[2][2] = { { 243., 274 }, { 243., 274 } };

		response[0] = normfactor[0] * parm[0][0] * exp(-y / parm[0][1]); // Right/bottom
		response[0] = 0.; // Only TOP readout
		response[1] = normfactor[1] * parm[1][0] * exp(-(pd_len - y) / parm[1][1]); // Left/top
		response[0] = response[1]; // Forcing ADC1 to have the only readout channel
		response[1] = 0.;

		//cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
		// cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
	}

//  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
// cout <<  " res[0]: " << response[0] << " channel " << channel<<endl;
	return response;
}

double veto_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks) {
//Example of Birk attenuation law in organic scintillators.
//adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
//
// Taken from GEANT4 examples advanced/amsEcal and extended/electromagnetic/TestEm3
//
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl);
	}
	return response;
}

double veto_HitProcess::BirksAttenuation2(double destep, double stepl, int charge, double birks) {
//Extension of Birk attenuation law proposed by Chou
// see G.V. O'Rielly et al. Nucl. Instr and Meth A368(1996)745
//
//
	double C = 9.59 * 1E-4 * mm * mm / MeV / MeV;
	double response = destep;
	if (birks * destep * stepl * charge != 0.) {
		response = destep / (1. + birks * destep / stepl + C * pow(destep / stepl, 2.));
	}
	return response;
}

map<string, vector<int> > veto_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map<string, vector<int> > MH;

	return MH;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> veto_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;

// first, identify the cells that would have electronic noise
// then instantiate hit with energy E, time T, identifier IDF:
//
// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

// push to noiseHits collection:
// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}

// - charge: returns charge/time digitized information / step
map<int, vector<double> > veto_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map<int, vector<double> > CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double veto_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

