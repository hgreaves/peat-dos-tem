#include "WildFire.h"

WildFire::WildFire() {
	VSMburn = 0.4; // a threshold value of VWC for burn

	r_ag2tot_cn = 0.8; //above-ground C to total
	r_bg2tot_cn = (1 - r_ag2tot_cn); //below-ground C to total

	r_live2ag_cn = 0.01; //live C of ag C
	r_dead2ag_cn = 0.76; //dead C of ag C
	r_burn2ag_cn = 1 - r_live2ag_cn - r_dead2ag_cn;

	r_retain_n = 0.3; //calculated from Harden et al., 2003 (ATHarden42003a)

	r_live2bg_cn = 0.01;
	//r_burn2bg_c;//determined every fire event
	// r_dead2bg_c;

}
;

WildFire::~WildFire() {

}
;

void WildFire::initializeParameter(const int & drgtypep, const int & vegtypep) {
	vegtype = vegtypep;
	drgtype = drgtypep;

	firpar.burnthick_max = chtlu->burnthick_max[vegtype];
	firpar.vegcombust = chtlu->vegcombust[vegtype];
	firpar.vegslash = chtlu->vegslash[vegtype];

}
;

void WildFire::initializeState() {
	fd->ysf = 0;
	fd->y_a2soi.orgn = 0;
}
;

void WildFire::initializeState5restart(RestartData *resin) {

	fd->ysf = resin->ysf; //resin->getYSF(fd->ysf, fd->cd->reschtid);
	fd->y_a2soi.orgn = resin->burnedn; //resin->getBURNEDN(fd->y_a2soi.orgn, fd->cd->reschtid);

}
;

//Yuan: modifying the following method, return the first fire year, if any
void WildFire::prepareDrivingData(const bool runsp, const bool runtr) {
	//initialize with -1
	for (int in = 0; in < MAX_FIR_OCR_NUM; in++) {
		occur[in] = -1;
		season[in] = -1;
		month[in] = -1;
		vegseverity[in] = -1;
		soiseverity[in] = -1;
		size[in] = -1;
	}

	//Yuan: season's month index order (0~11):
	int morder[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0 }; //Yuan: season: 0, 1(early fire), 2(late fire), and 3 with 3 months in the order
	vector<int> firemonths;

	int calyr = 0;

	firstfireyr = -1; //Yuan: first fire year specified in sp/tr: fire.nc

	if (runsp || runtr) {
		firstfireyr = END_TR_YR; // the latest possible year to have a fire

		//from sp-fire.nc
		for (int in = 0; in < MAX_SP_FIR_OCR_NUM; in++) {
			calyr = fd->cd->spfireyear[in];
			if (calyr != -1) { //Yuan: fire year may be BC, But -1 reserved for none

				if (firstfireyr >= calyr)
					firstfireyr = calyr;

				occur[in] = calyr - BEG_SP_YR; //Yuan: 1001 -> spin-up begining year
				season[in] = fd->cd->spseason[in];
				vegseverity[in] = fd->cd->spvegseverity[in];
				soiseverity[in] = fd->cd->spsoiseverity[in];

				int fsindx = season[in];
				for (int i = 0; i < 3; i++)
					firemonths.push_back(morder[fsindx * 3 + i]);
				random_shuffle(firemonths.begin(), firemonths.end()); //randomly pick-up a month for fire occurence
				month[in] = *firemonths.begin();
				firemonths.clear();

				//Yuan: the following modified based on the change of grid-level input:
				if (calyr < fd->gd->fireyear[0]) {
					size[in] = fd->gd->firesize[0]; //Yuan: need further modification here??
				} else {
					for (int ifs = 0; ifs < MAX_FSIZE_DRV_YR; ifs++) {
						if (calyr < fd->gd->fireyear[ifs]) {
							size[in] = fd->gd->firesize[ifs];
							break;
						}
					}
				}

			}
		}
	}

	if (runtr) {
		//from tr-fire.nc
		for (int in = 0; in < MAX_TR_FIR_OCR_NUM; in++) {
			calyr = fd->cd->trfireyear[in];
			if (calyr != -1) { //Yuan: fire year may be BC, But -1 reserved for none

				if (firstfireyr >= calyr)
					firstfireyr = calyr;

				occur[in + MAX_SP_FIR_OCR_NUM] = calyr - BEG_TR_YR; //Yuan: 1001 -> spin-up begining year
				season[in + MAX_SP_FIR_OCR_NUM] = fd->cd->trseason[in];
				vegseverity[in + MAX_SP_FIR_OCR_NUM]
						= fd->cd->trvegseverity[in];
				soiseverity[in + MAX_SP_FIR_OCR_NUM]
						= fd->cd->trsoiseverity[in];

				int fsindx = season[in + MAX_SP_FIR_OCR_NUM];
				for (int i = 0; i < 3; i++)
					firemonths.push_back(morder[fsindx * 3 + i]);
				random_shuffle(firemonths.begin(), firemonths.end()); //randomly pick-up a month for fire occurence
				month[in + MAX_SP_FIR_OCR_NUM] = *firemonths.begin();
				firemonths.clear();

				if (calyr < fd->gd->fireyear[0]) {
					size[in + MAX_SP_FIR_OCR_NUM] = fd->gd->firesize[0]; //Yuan: need further modification here??
				} else {
					for (int ifs = 0; ifs < MAX_FSIZE_DRV_YR; ifs++) {
						if (calyr < fd->gd->fireyear[ifs]) {
							size[in + MAX_SP_FIR_OCR_NUM]
									= fd->gd->firesize[ifs];
							break;
						}
					}
				}

			}
		}
	}

}
;

//Yuan: the fire occurrence (and data) is input (cohort-level info),or FRI derived (grid-level info)
//Yuan: almost rewriting the code
int WildFire::getOccur(const int &yrind, const int & mind,
		const bool & friderived) {
	int occ = 0;

	if (friderived) {
		if (yrind % fd->gd->fri == 0 && yrind > 0) {
			if (mind == (fd->gd->fireseason[0]) * 3 + 2) { //Yuan: the middle month in gd->fireseason
				occ = 1;
				return occ;
			}
		}

	} else {
		for (int i = 0; i < MAX_FIR_OCR_NUM; i++) {
			if (occur[i] == yrind) {
				if (mind == month[i]) {
					occ = 1;
					return occ;
				}
			}
		}
	}

	return occ;

}
;

//Yuan: the fire occurrence (and data) is input or FRI derived
void WildFire::burn(const int & yrind, const bool & friderived) {
	fd->burn();

	if (!friderived) {
		//Yuan: the following modified for less memory consumption
		for (int i = 0; i < MAX_FIR_OCR_NUM; i++) { // Yuan: Fire data reads from 4 arrays (cohort-level info)
			if (occur[i] == yrind) {
				onesize = size[i];
				oneseason = season[i];
				onevegseverity = vegseverity[i];
				onesoiseverity = soiseverity[i];

				break;
			}

		}

	} else {
		if (yrind % fd->gd->fri == 0) { //Yuan: occurrence based on FRI (grid-level info)
			onesize = fd->gd->firesize[0]; //temporarily set
			oneseason = fd->gd->fireseason[0];
			onevegseverity = -1;
			onesoiseverity = -1;
		}

	}

	// for soil part // THERMOKARST VERSION
	updateBurnThickness();
	double burndepth = fd->y_soid.burnthick;
	cout << "Wildfire::burn -- burndepth at start of burn: " << burndepth << "\n";
	double totbotdepth = 0.;
	double burnedsolc = 0.;
	r_burn2bg_cn = 0.; //initialize

	for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {

		totbotdepth += ed->m_sois.dz[il];

		if (totbotdepth <= burndepth) { //remove all the orgc in this layer, consider nitrogen later
			burnedsolc += bd->m_sois.reac[il] + bd->m_sois.nonc[il];
			bd->m_sois.reac[il] = 0.;
			bd->m_sois.nonc[il] = 0.;

			r_burn2bg_cn += ed->m_sois.rootfrac[il];
			ed->m_sois.rootfrac[il] = 0.;
		} else {
			double partleft = totbotdepth - burndepth;
			//calculate the left carbon
			if (partleft < ed->m_sois.dz[il]) {
				if (ed->m_sois.type[il] <= 2) { 
					burnedsolc += (1 - partleft / ed->m_sois.dz[il])
							* (bd->m_sois.reac[il] + bd->m_sois.nonc[il]);
					bd->m_sois.reac[il] *= partleft / ed->m_sois.dz[il];
					bd->m_sois.nonc[il] *= partleft / ed->m_sois.dz[il];

				r_burn2bg_cn += (1 - partleft / ed->m_sois.dz[il])
						* ed->m_sois.rootfrac[il];
				ed->m_sois.rootfrac[il] *= partleft / ed->m_sois.dz[il];

				} else {
				//nothing can happen here
					break;
				}
		

			} else {
				break;
			}

		}
	}


// Heather's note: the following preliminary code chunk was based on the original code (following), but was 
// subsequently split (above) to put the thickness calculation (bthick) into getBurnThick
// and keep this part as just the carbon calculation,  to better mimic the original code structure

	// // Modified code for thermokarst carbon/layer loss (copied from below and modified)
	// double solz = 0.0; // overall org layer thickness
	// double solztk = 0.0; // org layer thickness lost to thermokarst 
	// double solc = 0.0;  // overall org layer carbon 
	// double cumc = 0.0;  // cumulative org layer carbon (while progressing through loop)

	// for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {
	// 	if (ed->m_sois.type[il] <= 2) { //moss is zero, shlw peat is 1 and deep org is 2
	// 			solc += bd->m_sois.reac[il] + bd->m_sois.nonc[il]; // retrieve total org soil carbon
	// 			solz += ed->m_sois.dz[il]; // retrieve total org soil depth
	// 	}
	// }
	// cout << "Total org layer carbon (solc before tk): " << solc << "\n";

	// double solctk = solc * 0.30; // org layer carbon to be lost to thermokarst
	
	// cout << "Org layer carbon to be lost to tk (solctk): " << solctk << "\n";
	// cout << "Total org layer depth (solz): " << solz << "\n";

	// for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {
	// 	if (ed->m_sois.type[il] <= 2) { //moss is zero, shlw peat is 1 and deep org is 2
	// 		// retrieve and remove carbon to be lost from org layer soils to thermokarst
	// 		cumc += bd->m_sois.reac[il] + bd->m_sois.nonc[il]; //accumulate soil carbon
	// 		cout <<"cumc: " << cumc << "\n";
	// 		if (cumc < solctk){
	// 			burnedsolc += bd->m_sois.reac[il] + bd->m_sois.nonc[il]; //add layer carbon to lost carbon
	// 			solztk += ed->m_sois.dz[il]; // add layer depth to lost depth
				
	// 			cout << "cumc < solc: burnedsolc: " << burnedsolc << "\n";
	// 			cout << "solztk: " << solztk << "\n";
				
	// 			bd->m_sois.reac[il] = 0.; // set layer carbon to zero
	// 			bd->m_sois.nonc[il] = 0.;

	// 			cout << "r_burn2bg_cn: " << r_burn2bg_cn << "\n"; 
	// 			r_burn2bg_cn += ed->m_sois.rootfrac[il]; // add layer roots to lost carbon (?) and set to zero
	// 			ed->m_sois.rootfrac[il] = 0.;
	// 		} else {
	// 			double cleft = cumc - solctk; // carbon remaining in layer
	// 			cout << "Partial layer reached; cleft: "<< cleft << "\n";
	// 			ed->m_sois.dz[il] *= cleft / (bd->m_sois.reac[il] + bd->m_sois.nonc[il]); // reduce thickness proportionately
	// 			solztk += ed->m_sois.dz[il]; // add partial thickness to lost depth
				
	// 			if (cleft < (bd->m_sois.reac[il] + bd->m_sois.nonc[il])) {
	// 				if ((solz - solztk) < 0.02){
	// 					cleft *= ((0.02 - (solz - solztk))/ed->m_sois.dz[il]);
	// 					ed->m_sois.dz[il] = 0.02 - (solz - solztk);
	// 				}
	// 				cout << "dz: " << ed->m_sois.dz[il] <<"\n";
	// 				cout << "cleft: " << cleft <<"\n";
	// 				//calculate the left carbon
	// 				burnedsolc += (bd->m_sois.reac[il] + bd->m_sois.nonc[il]) - cleft ;
	// 				cout << "cumc >= solc: burnedsolc: " << burnedsolc << "\n";
	// 				cout << "solztk: " << solztk << "\n";
	// 				bd->m_sois.reac[il] *= cleft / (bd->m_sois.reac[il] + bd->m_sois.nonc[il]);//0;//lwcumorgc - upcumorgc;
	// 				bd->m_sois.nonc[il] *= cleft / (bd->m_sois.reac[il] + bd->m_sois.nonc[il]);//0;//lwcumorgc - upcumorgc;
						
	// 				r_burn2bg_cn += (1 - (cleft / (bd->m_sois.reac[il] + bd->m_sois.nonc[il])))
	// 							* ed->m_sois.rootfrac[il];
	// 				ed->m_sois.rootfrac[il] *= (cleft / (bd->m_sois.reac[il] + bd->m_sois.nonc[il]));

				
	// 			} else {
	// 				break;
	// 			}

	// 		}
			
	// 	} else {
	// 		break;
	// 	}
	// }


/// ORIGINAL VERSION:

	// for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {

	// 		totbotdepth += ed->m_sois.dz[il];

	// 		if (totbotdepth <= burndepth) { //remove all the orgc in this layer, consider nitrogen later
	// 			burnedsolc += bd->m_sois.reac[il] + bd->m_sois.nonc[il];
	// 			bd->m_sois.reac[il] = 0.;
	// 			bd->m_sois.nonc[il] = 0.;

	// 			r_burn2bg_cn += ed->m_sois.rootfrac[il];
	// 			ed->m_sois.rootfrac[il] = 0.;
	// 		} else {
	// 			double partleft = totbotdepth - burndepth;
	// 			//calculate the left carbon
	// 			if (partleft < ed->m_sois.dz[il]) {
	// 				if (ed->m_sois.type[il] == 1) { //shallow organic
	// 					//double upcumorgc = getCumulativeCarbonBD(burndepth);
	// 					//double lwcumorgc = getCumulativeCarbonBD(totbotdepth);
	// 					burnedsolc += (1 - partleft / ed->m_sois.dz[il])
	// 							* (bd->m_sois.reac[il] + bd->m_sois.nonc[il]);
	// 					bd->m_sois.reac[il] *= partleft / ed->m_sois.dz[il];//0;//lwcumorgc - upcumorgc;
	// 					bd->m_sois.nonc[il] *= partleft / ed->m_sois.dz[il];//0;//lwcumorgc - upcumorgc;

	// 				} else { //deep organic
	// 					burnedsolc += (1 - partleft / ed->m_sois.dz[il])
	// 							* (bd->m_sois.reac[il] + bd->m_sois.nonc[il]);
	// 					bd->m_sois.reac[il] *= partleft / ed->m_sois.dz[il];
	// 					bd->m_sois.nonc[il] *= partleft / ed->m_sois.dz[il];

	// 				}

	// 				r_burn2bg_cn += (1 - partleft / ed->m_sois.dz[il])
	// 						* ed->m_sois.rootfrac[il];
	// 				ed->m_sois.rootfrac[il] *= partleft / ed->m_sois.dz[il];

	// 			} else {
	// 				//nothing can happen here
	// 				break;
	// 			}
			

	// 		} else {
	// 			break;
	// 		}

	// }

	//soil N burned and returned	
	double burnedsoln = burnedsolc / chtlu->cnsoil[vegtype];

	double vola_soln = burnedsoln * (1 - r_retain_n);
	double reta_soln = burnedsoln * r_retain_n;

	// for vegetation
	// above ground

	if (fd->useseverity) {
		r_burn2ag_cn = onevegseverity;
		r_dead2ag_cn = 1 - r_burn2ag_cn - r_live2ag_cn;
	}

	double comb_ag_vegc = bd->m_vegs.c * r_ag2tot_cn * r_burn2ag_cn;
	double live_ag_vegc = bd->m_vegs.c * r_ag2tot_cn * r_live2ag_cn;
	double dead_ag_vegc = bd->m_vegs.c * r_ag2tot_cn * r_dead2ag_cn;

	double comb_ag_strn = bd->m_vegs.strn * r_ag2tot_cn * r_burn2ag_cn;
	double live_ag_strn = bd->m_vegs.strn * r_ag2tot_cn * r_live2ag_cn;
	double dead_ag_strn = bd->m_vegs.strn * r_ag2tot_cn * r_dead2ag_cn;
	double vola_ag_strn = comb_ag_strn * (1 - r_retain_n);
	double reta_ag_strn = comb_ag_strn * (r_retain_n);

	double comb_ag_ston = bd->m_vegs.ston * r_ag2tot_cn * r_burn2ag_cn;
	double live_ag_ston = bd->m_vegs.ston * r_ag2tot_cn * r_live2ag_cn;
	double dead_ag_ston = bd->m_vegs.ston * r_ag2tot_cn * r_dead2ag_cn;
	double vola_ag_ston = comb_ag_ston * (1 - r_retain_n);
	double reta_ag_ston = comb_ag_ston * (r_retain_n);

	// below ground, the burn ratio is calculated at each fire event above
	if (r_burn2bg_cn > 0.99)
		r_burn2bg_cn = 0.99;//there are at lest 0.01 live bg vegc
	r_dead2bg_cn = 1 - r_burn2bg_cn - r_live2bg_cn;
	double comb_bg_vegc = bd->m_vegs.c * r_bg2tot_cn * r_burn2bg_cn;
	double live_bg_vegc = bd->m_vegs.c * r_bg2tot_cn * r_live2bg_cn;
	double dead_bg_vegc = bd->m_vegs.c * r_bg2tot_cn * r_dead2bg_cn;

	double comb_bg_strn = bd->m_vegs.strn * r_bg2tot_cn * r_burn2bg_cn;
	double live_bg_strn = bd->m_vegs.strn * r_bg2tot_cn * r_live2bg_cn;
	double dead_bg_strn = bd->m_vegs.strn * r_bg2tot_cn * r_dead2bg_cn;
	double vola_bg_strn = comb_bg_strn * (1 - r_retain_n);
	double reta_bg_strn = comb_bg_strn * (r_retain_n);

	double comb_bg_ston = bd->m_vegs.ston * r_bg2tot_cn * r_burn2bg_cn;
	double live_bg_ston = bd->m_vegs.ston * r_bg2tot_cn * r_live2bg_cn;
	double dead_bg_ston = bd->m_vegs.ston * r_bg2tot_cn * r_dead2bg_cn;
	double vola_bg_ston = comb_bg_ston * (1 - r_retain_n);
	double reta_bg_ston = comb_bg_ston * (r_retain_n);

	/// output the biome-atm exchange data during Fire and BGC
	// for soil organic layers
	fd->y_soi2a.orgc = burnedsolc;
	cout << "burnedsolc 377: " << burnedsolc << "\n";
	fd->y_soi2a.orgn = vola_soln; //nitrogen emission due to burn

	if (fd->y_soi2a.orgn > bd->m_sois.orgn) {
		fd->y_soi2a.orgn = bd->m_sois.orgn;
	}

	// then for vegetation

	fd->y_v2a.orgc = comb_ag_vegc + comb_bg_vegc;
	fd->y_v2a.orgn = vola_bg_ston + vola_bg_strn + vola_ag_ston + vola_ag_strn;

	bd->m_vegs.deadc = dead_ag_vegc;
	bd->m_vegs.deadn = dead_bg_strn + dead_ag_ston + dead_ag_strn
			+ dead_bg_ston;
	bd->m_vegs.c = live_bg_vegc + live_ag_vegc;
	if (bd->m_vegs.c < 1)
		bd->m_vegs.c = 1;
	bd->m_vegs.strn = live_bg_strn + live_ag_strn;
	bd->m_vegs.ston = live_bg_ston + live_ag_ston;

	// calculate the amount of nitrogen returned into soil system every year

	bd->m_sois.orgn += (reta_bg_ston + reta_bg_strn + reta_ag_ston
			+ reta_ag_strn + reta_soln - vola_soln);
	//need to input the nitrogen from wood debris
	bd->m_sois.orgn += bd->m_vegs.deadn;

	fd->y_a2soi.orgn = (fd->y_soi2a.orgn + fd->y_v2a.orgn) / fd->gd->fri; //Yuan: ???

	//put the dead bg C into each layer

	for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {
		if (ed->m_sois.rootfrac[il] >= 0.001) {
			bd->m_sois.nonc[il] += dead_bg_vegc * ed->m_sois.rootfrac[il];
		}
	}

	//also reset bd.foliagemx
	bd->foliagemx = 0;

}
;

void WildFire::updateBurnThickness() {
	//for initial test, assume only 10 cm will be left
	double bdepth = 0.;
	//get the total orgnic depth
	//////////////////////////////////
	///Rule 1: only organic layer can be burned (Site Related)
	///Rule 2: should be less than the burn thickness max (Cohort Related)
	///Rule 3: should be less than the water table depth

	// Heather says: remove constraint to test thermokarst, replace with standard calculation of totorg
	// double totorg = 0.;

	// for (int i = 0; i < ed->m_soid.actual_num_soil; i++) {
		
	// 	if (ed->m_sois.type[i] <= 2 && ed->m_soid.allvwc[i] <= VSMburn) { //Yuan: VSM constrained burn thickness
	// 		totorg += ed->m_sois.dz[i];
	// 	} else {
	// 		break;
	// 	}

	// }

	double totorg = ed->m_soid.mossthick + ed->m_soid.shlwthick
			+ ed->m_soid.deepthick;
	
	bdepth = getBurnThick();
	if (bdepth > totorg) {
		bdepth = totorg;
	} //Yuan:

	fd->y_soid.burnthick = bdepth;
	cout << "bdepth (bthick) at the end of updateBurnThickness: " << bdepth << "\n";
}
;

double WildFire::getBurnThick() {
	cout << "getBurnThick: \n";
	double bthick = 0; // this is currently the thickness lost to thermokarst
	double rtk = 0.30; // fraction of soil c loss to thermokarst
	// double maxburn = firpar.burnthick_max;
	double totorg = ed->m_soid.mossthick + ed->m_soid.shlwthick
			+ ed->m_soid.deepthick;

	double soldz = 0.0; // overall org layer thickness
	double solc = 0.0;  // overall org layer carbon 
	double cumc = 0.0;  // cumulative org layer carbon (while progressing through loop)
	double solctk = 0.0; // how much org layer c is lost (just to check against target c loss)

	// Get total sol carbon and thickness (we already have thickness with totorg... whatever.)
	for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {
		if (ed->m_sois.type[il] <= 2) { //moss is zero, shlw peat is 1 and deep org is 2
				solc += bd->m_sois.reac[il] + bd->m_sois.nonc[il]; // retrieve total org soil carbon
				soldz += ed->m_sois.dz[il]; // retrieve total org soil depth
		}
	}
	cout << "Total org layer thickness (totorg): " << totorg << "\n";
	cout << "Total org layer thickness (soldz): " << soldz << "\n";
	cout << "Total org layer carbon (solc before tk): " << solc << "\n";
	
	double solctk_target = solc * rtk; // org layer carbon to be lost to thermokarst
	
	cout << "Total c loss target (solctk_target): " << solctk_target << "\n";

	for (int il = 0; il < ed->m_soid.actual_num_soil; il++) {
		if (ed->m_sois.type[il] <= 2) { //moss is zero, shlw peat is 1 and deep org is 2
			// retrieve carbon to be lost from org layer soils to thermokarst and convert to burn-like thickness
			cumc += bd->m_sois.reac[il] + bd->m_sois.nonc[il]; //accumulate soil carbon
			cout << "Layer: " << il <<" Accumulated soil c so far (cumc): " << cumc << "\n";
			if (cumc < solctk_target){
				bthick += ed->m_sois.dz[il]; // add layer depth to lost depth
				cout << "bthick: " << bthick << "\n";
				
			} else {
				double cleft = cumc - solctk_target; // carbon that will remain in layer
				cout << "Partial layer reached; cleft: " << cleft << "\n";
					
				if (cleft < (bd->m_sois.reac[il] + bd->m_sois.nonc[il])){
					double dzlost = ed->m_sois.dz[il] * (1 - (cleft / (bd->m_sois.reac[il] + bd->m_sois.nonc[il]))); 
					cout << "dz of partial layer: " << ed->m_sois.dz[il] << "\n";
					cout << "dzlost from partial layer: " << dzlost << "\n";
					bthick += dzlost; // add partial thickness to lost depth
					cout << "Total bthick: " << bthick << "\n";

					solctk = cumc - (bd->m_sois.reac[il] + bd->m_sois.nonc[il]) + (bd->m_sois.reac[il] + bd->m_sois.nonc[il]) * (dzlost / ed->m_sois.dz[il]);
				} else {
					break;
				}
			}
			
		} else {
			break;
		}
	}

	cout << "getBurnThick Check: Target C loss: " << solctk_target << ", Actual C loss: " << solctk << "\n";
	
	// Commented out for thermokarst
	// // currently use half of the total organic
	// if (!fd->useseverity) {
	// 	if (fd->cd->drgtype == 0) {
	// 		if (oneseason == 0 || oneseason == 1 || oneseason == 3) { //Yuan: bug here? (fireseason: 0, 1, 2, 3)
	// 		//     			if(onesize<1){                     //Yuan: bug here? (firesize: 0, 1, 2, 3)
	// 			if (onesize <= 1) {
	// 				bthick = 0.54 * totorg;
	// 			} else {
	// 				bthick = 0.69 * totorg;
	// 			}
	// 		} else if (oneseason == 2) { //late season
	// 			bthick = 0.8 * totorg;
	// 		}
	// 		rtk = 0.0;
	// 	} else if (fd->cd->drgtype == 1) {
	// //		bthick = 0.48 * totorg;
	// 		rtk = 0.30;
	// 	}

	// } else {
	// 	bthick = onesoiseverity * totorg;
	// }

	cout << "bthick @ 593 after modified loop 593: " << bthick << "\n";
	if (bthick <= 0) {
		cout << "burn thickness: " << bthick << "\n";
		string msg = "burn thickness should be greater than zero";
		char* msgc = const_cast<char*> (msg.c_str());
		throw Exception(msgc, I_BURN_ZERO);
	}
	cout << "bthick before moss check: " << bthick << "\n";
	cout << "mossthick: " << ed->m_soid.mossthick << "\n";
	if (bthick < ed->m_soid.mossthick) {
		bthick = ed->m_soid.mossthick;
	}
 	cout << "bthick after moss check: " << bthick << "\n";
 	cout << "totorg @ 605: " << totorg << "\n";
	if (totorg - bthick < 0.02) { //there are at least 2 cm humic orgnanic left
		bthick = totorg - 0.02;
	}
	cout << "bthick at end of getBurnThick: " << bthick << "\n";
	return bthick;
}
;

void WildFire::setCohortLookup(CohortLookup* chtlup) {
	chtlu = chtlup;
}
;

void WildFire::setEnvData(EnvData* edp) {
	ed = edp;
}
;

void WildFire::setBgcData(BgcData* bdp) {
	bd = bdp;
}
;

void WildFire::setFirData(FirData* fdp) {
	fd = fdp;
}

