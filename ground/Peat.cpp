#include "Peat.h"

Peat::Peat() {

}

Peat::~Peat() {

}

bool Peat::sameShlwThickRange(const double & thickness) {

	/// if the thickness and shlwthick not in the same range, there is a need to reconstruct peat layers
	bool same = false;
	if (shlwnum == 2) {
		if (shlwdza[0] > shlwdza[1]) {
			same = false;
			return same;
		}

	} else if (shlwnum == 3) {
		if (shlwdza[0] > shlwdza[1] || shlwdza[0] > shlwdza[2] || shlwdza[1]
				> shlwdza[2]) {
			same = false;
			return same;
		}
	}

	if (thickness < 0.02 && shlwthick < 0.02) {
		return true;
	} else if (thickness < 0.04 && thickness >= 0.01 && shlwthick < 0.04
			&& shlwthick >= 0.02) {
		return true;
	} else if (thickness < 0.06 && thickness >= 0.04 && shlwthick < 0.06
			&& shlwthick >= 0.04) {
		return true;
	} else if (thickness < 0.10 && thickness >= 0.06 && shlwthick < 0.10
			&& shlwthick >= 0.06) {
		return true;
	} else if (thickness < 0.15 && thickness >= 0.10 && shlwthick < 0.15
			&& shlwthick >= 0.10) {
		return true;
	} else if (thickness < 0.20 && thickness >= 0.15 && shlwthick < 0.20
			&& shlwthick >= 0.15) {
		return true;

	} else if (thickness < 0.28 && thickness >= 0.20 && deepthick < 0.28
			&& deepthick >= 0.20) {
		return true;
	} else if (thickness < 0.40 && thickness >= 0.28 && deepthick < 0.40
			&& deepthick >= 0.28) {
		return true;
	} else if (thickness >= 0.4 && deepthick >= 0.4) {
		return true;
	}

	return same;
}

bool Peat::sameDeepThickRange(const double & thickness) {
	/// if the thickness and deepthick not in the same range, there is a need to reconstruct peat layers
	bool same = false;

	if (shlwchanged) {
		same = false;
		return same;
	}

	int newdeepnum = 0;

	if (thickness <= dlstshlw)
		newdeepnum = 1;
	else if (thickness > dlstshlw && thickness <= 2 * dlstshlw)
		newdeepnum = 2;
	else if (thickness > 2 * dlstshlw && thickness <= 3 * dlstshlw)
		newdeepnum = 3;
	else if (thickness > 3 * dlstshlw && thickness <= 4 * dlstshlw)
		newdeepnum = 4;
	else if (thickness > 4 * dlstshlw && thickness <= 5 * dlstshlw)
		newdeepnum = 5;
	else if (thickness > 5 * dlstshlw && thickness <= 6 * dlstshlw)
		newdeepnum = 6;
	else if (thickness > 6 * dlstshlw && thickness <= 7 * dlstshlw)
		newdeepnum = 7;
	else if (thickness > 7 * dlstshlw && thickness <= 8 * dlstshlw)
		newdeepnum = 8;
	else if (thickness > 8 * dlstshlw && thickness <= 9 * dlstshlw)
		newdeepnum = 9;
	else
		newdeepnum = 10;

	/*
	 if (thickness < 3 * dlstshlw) {

	 newdeepnum = 1;
	 } else if (thickness >= 3 * dlstshlw && thickness < 6 * dlstshlw) {

	 newdeepnum = 2;
	 } else {

	 newdeepnum = 3;

	 }
	 */
	if (newdeepnum == deepnum) {
		same = true;
	}

	return same;
}

void Peat::initShlwThicknesses(const double & thickness) {

	oldshlwthick = shlwthick;
	if (shlwnum > 0)
		olddlstshlw = shlwdza[shlwnum - 1]; //Yuan:???
	shlwthick = thickness;

	if (thickness < 0.02)
		shlwthick = 0.02;
	shlwnum = 0;
	if (thickness < 0.00) { //Yuan: always allow at least ONE shlwnum
		shlwdza[0] = -1;
		shlwdza[1] = -1;
		shlwdza[2] = -1;
		shlwnum = 0;
	} else if (thickness < 0.04) {
		shlwdza[0] = thickness;
		shlwdza[1] = -1;
		shlwdza[2] = -1;
		shlwnum = 1;
	} else if (thickness < 0.06) {
		shlwdza[0] = 0.02;
		shlwdza[1] = thickness - 0.02;
		shlwdza[2] = -1;
		shlwnum = 2;
	} else if (thickness < 0.10) {
		shlwdza[0] = 0.03;
		shlwdza[1] = thickness - 0.03;
		shlwdza[2] = -1;
		shlwnum = 2;
	} else if (thickness < 0.15) {
		shlwdza[0] = 0.02;
		shlwdza[1] = 0.04;
		shlwdza[2] = thickness - 0.06;
		shlwnum = 3;
	} else if (thickness < 0.20) {
		shlwdza[0] = 0.03;
		shlwdza[1] = 0.06;
		shlwdza[2] = thickness - 0.09;
		shlwnum = 3;
	} else if (thickness < 0.28) {
		shlwdza[0] = 0.04;
		shlwdza[1] = 0.08;
		shlwdza[2] = thickness - 0.12;
		shlwnum = 3;
	} else if (thickness < 0.40) {
		shlwdza[0] = 0.05;
		shlwdza[1] = 0.11;
		shlwdza[2] = thickness - 0.16;
		shlwnum = 3;
	} else if (thickness >= 0.4) {

		shlwdza[0] = 0.1;
		shlwdza[1] = 0.2;
		shlwdza[2] = thickness - 0.3;
		shlwnum = 3;
	} else {
		//cout << "the thickness of shlw organic is too big " << thickness<< "\n";
		//exit(-1);
	}
	shlwthick = thickness;
	dlstshlw = shlwdza[shlwnum - 1];
}

void Peat::initDeepThicknesses(const double & thickness) {
	olddeepthick = deepthick;

	deepthick = thickness;
	if (thickness < 0.0)
		deepthick = 0.0;

	deepnum = 0;

	if (dlstshlw > 0) {
		if (deepthick <= dlstshlw) {
			deepdza[0] = deepthick;
			deepdza[1] = -1;
			deepdza[2] = -1;
			deepdza[3] = -1;
			deepdza[4] = -1;
			deepdza[5] = -1;
			deepdza[6] = -1;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 1;
		} else if (deepthick > dlstshlw && deepthick <= 2 * dlstshlw) {
			deepdza[0] = 1. / 3. * deepthick;
			deepdza[1] = 2. / 3. * deepthick;
			deepdza[2] = -1;
			deepdza[3] = -1;
			deepdza[4] = -1;
			deepdza[5] = -1;
			deepdza[6] = -1;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 2;
		} else if (deepthick > 2 * dlstshlw && deepthick <= 3 * dlstshlw) {
			deepdza[0] = 1. / 6. * deepthick;
			deepdza[1] = 2. / 6. * deepthick;
			deepdza[2] = 3. / 6. * deepthick;
			deepdza[3] = -1;
			deepdza[4] = -1;
			deepdza[5] = -1;
			deepdza[6] = -1;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 3;
		} else if (deepthick > 3 * dlstshlw && deepthick <= 4 * dlstshlw) {
			deepdza[0] = 1. / 10. * deepthick;
			deepdza[1] = 2. / 10. * deepthick;
			deepdza[2] = 3. / 10. * deepthick;
			deepdza[3] = 4. / 10. * deepthick;
			deepdza[4] = -1;
			deepdza[5] = -1;
			deepdza[6] = -1;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 4;
		} else if (deepthick > 4 * dlstshlw && deepthick <= 5 * dlstshlw) {
			deepdza[0] = 1. / 15. * deepthick;
			deepdza[1] = 2. / 15. * deepthick;
			deepdza[2] = 3. / 15. * deepthick;
			deepdza[3] = 4. / 15. * deepthick;
			deepdza[4] = 5. / 15. * deepthick;
			deepdza[5] = -1;
			deepdza[6] = -1;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 5;
		} else if (deepthick > 5 * dlstshlw && deepthick <= 6 * dlstshlw) {
			deepdza[0] = 1. / 21. * deepthick;
			deepdza[1] = 2. / 21. * deepthick;
			deepdza[2] = 3. / 21. * deepthick;
			deepdza[3] = 4. / 21. * deepthick;
			deepdza[4] = 5. / 21. * deepthick;
			deepdza[5] = 6. / 21. * deepthick;
			deepdza[6] = -1;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 6;
		} else if (deepthick > 6 * dlstshlw && deepthick <= 7 * dlstshlw) {
			deepdza[0] = 1. / 28. * deepthick;
			deepdza[1] = 2. / 28. * deepthick;
			deepdza[2] = 3. / 28. * deepthick;
			deepdza[3] = 4. / 28. * deepthick;
			deepdza[4] = 5. / 28. * deepthick;
			deepdza[5] = 6. / 28. * deepthick;
			deepdza[6] = 7. / 28. * deepthick;
			deepdza[7] = -1;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 7;
		} else if (deepthick > 7 * dlstshlw && deepthick <= 8 * dlstshlw) {
			deepdza[0] = 1. / 36. * deepthick;
			deepdza[1] = 2. / 36. * deepthick;
			deepdza[2] = 3. / 36. * deepthick;
			deepdza[3] = 4. / 36. * deepthick;
			deepdza[4] = 5. / 36. * deepthick;
			deepdza[5] = 6. / 36. * deepthick;
			deepdza[6] = 7. / 36. * deepthick;
			deepdza[7] = 8. / 36. * deepthick;
			deepdza[8] = -1;
			deepdza[9] = -1;

			deepnum = 8;
		} else if (deepthick > 8 * dlstshlw && deepthick <= 9 * dlstshlw) {
			deepdza[0] = 1. / 45. * deepthick;
			deepdza[1] = 2. / 45. * deepthick;
			deepdza[2] = 3. / 45. * deepthick;
			deepdza[3] = 4. / 45. * deepthick;
			deepdza[4] = 5. / 45. * deepthick;
			deepdza[5] = 6. / 45. * deepthick;
			deepdza[6] = 7. / 45. * deepthick;
			deepdza[7] = 8. / 45. * deepthick;
			deepdza[8] = 9. / 45. * deepthick;
			deepdza[9] = -1;

			deepnum = 9;
		} else {
			deepdza[0] = 1. / 55. * deepthick;
			deepdza[1] = 2. / 55. * deepthick;
			deepdza[2] = 3. / 55. * deepthick;
			deepdza[3] = 4. / 55. * deepthick;
			deepdza[4] = 5. / 55. * deepthick;
			deepdza[5] = 6. / 55. * deepthick;
			deepdza[6] = 7. / 55. * deepthick;
			deepdza[7] = 8. / 55. * deepthick;
			deepdza[8] = 9. / 55. * deepthick;
			deepdza[9] = 10. / 55. * deepthick;

			deepnum = 10;
		}
	} else {
		deepdza[0] = 1. / 55. * deepthick;
		deepdza[1] = 2. / 55. * deepthick;
		deepdza[2] = 3. / 55. * deepthick;
		deepdza[3] = 4. / 55. * deepthick;
		deepdza[4] = 5. / 55. * deepthick;
		deepdza[5] = 6. / 55. * deepthick;
		deepdza[6] = 7. / 55. * deepthick;
		deepdza[7] = 8. / 55. * deepthick;
		deepdza[8] = 9. / 55. * deepthick;
		deepdza[9] = 10. / 55. * deepthick;

		deepnum = 10;
		/*
		deepdza[0] = deepthick;
		deepdza[1] = -1;
		deepdza[2] = -1;
		deepdza[3] = -1;
		deepdza[4] = -1;
		deepdza[5] = -1;
		deepdza[6] = -1;
		deepdza[7] = -1;
		deepdza[8] = -1;
		deepdza[9] = -1;

		deepnum = 1;
		*/
	}
	/*
	 if (dlstshlw > 0) {
	 if (deepthick < 3 * dlstshlw) {
	 deepdza[0] = deepthick;
	 deepdza[1] = -1;
	 deepdza[2] = -1;
	 deepnum = 1;
	 } else if (deepthick >= 3 * dlstshlw && deepthick < 6 * dlstshlw) {
	 deepdza[0] = 1. / 3. * deepthick;
	 deepdza[1] = 2. / 3. * deepthick;
	 deepdza[2] = -1;
	 deepnum = 2;
	 } else {
	 deepdza[0] = 1. / 6. * deepthick;
	 deepdza[1] = 2. / 6. * deepthick;
	 deepdza[2] = 3. / 6. * deepthick;
	 deepnum = 3;
	 }
	 } else {
	 deepdza[0] = deepthick;
	 deepdza[1] = -1;
	 deepdza[2] = -1;
	 deepnum = 1;
	 }*/

}

void Peat::updateShlwThicknesses(int type[], double dz[], const int & maxnum) {
	shlwnum = 0;
	shlwthick = 0;
	for (int i = 0; i < maxnum; i++) {
		if (type[i] == 1) {
			shlwdza[shlwnum] = dz[i];
			shlwnum++;
			shlwthick += dz[i];
		} else {
			if (type[i] > 1) {
				break;
			}
		}
	}
}

void Peat::updateDeepThicknesses(int type[], double dz[], const int & maxnum) {
	deepnum = 0;
	deepthick = 0;

	for (int i = 0; i < maxnum; i++) {
		if (type[i] == 2) {
			deepdza[deepnum] = dz[i];
			deepnum++;
			deepthick += dz[i];
		} else {
			if (type[i] > 2) {
				break;
			}
		}
	}
	if (deepnum == 1) {
		if (deepthick < 0.02) {
			deepthick = 0.02;
			deepdza[0] = 0.02;
		}
	} else if (deepnum == 0) {
		deepnum = 1;
		deepthick = 0.02;
		deepdza[0] = 0.02;

	}

}
