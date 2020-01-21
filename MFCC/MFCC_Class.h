#pragma once

//#ifndef MFCC_H_
//#define MFCC_H_

#define _USE_MATH_DEFINES
#include<complex>
#include<vector>
#include<cmath>

#include "audioConstants.h"

typedef std::complex<double> complexD;

class MFCC_Class
{
	MFCC_Class();
	static std::vector<complexD> FFT(std::vector<complexD> frame);
	static complexD getWi(int i, int n);
	static std::vector<double> windowHamming(std::vector<complexD> frame);

	static std::vector<std::vector<double> > getMelFilters(int mfccCount, int frameLength, int frequency, int frequencyMin, int frequencyMax);

	static std::vector<double> logEnergySpectrum(std::vector<double> spectrum, std::vector<std::vector<double> > melFilters, int mfccCount);

	static std::vector<double> DCT(std::vector<double> data);

	static double convertToMel(double f)
	{
		return 1127.0 * std::log(1.0 + f / 700.0);
	}
	static double convertFromMel(double m)
	{
		return 700.0 * (std::exp(m / 1127.0) - 1.0);
	}
	static double superpose(double f, int frameLength, int frequency)
	{
		return std::floor((frameLength + 1) * f / (double)frequency);
	}

public:
	static std::vector<double> computeMFCC(std::vector<double> frame, int frequency, int mfccCount = MFCC_COUNT, int frequencyMin = FRECUENCY_MIN, int frequencyMax = FRECUENCY_MAX);
};

//#endif /* MFCC_H_ */


