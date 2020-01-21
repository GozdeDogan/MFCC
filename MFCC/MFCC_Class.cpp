#include "MFCC_Class.h"

std::vector<complexD> MFCC_Class::FFT(std::vector<complexD> frame)
{
	int n = frame.size();
	if (n == 1)
	{
		return frame;
	}

	std::vector<complexD> x0(n / 2), x1(n / 2);
	std::vector<complexD> w(n);
	double alpha;

	for (int i = 0, j = 0; i < n; i += 2, j++)
	{
		x0[j] = frame[i];
		x1[j] = frame[i + 1];

		w[i] = getWi(i, n);
		w[i + 1] = getWi(i + 1, n);
	}

	std::vector<complexD> y0 = FFT(x0), y1 = FFT(x1), y(n);

	for (int i = 0; i < n; i++)
	{
		y[i] = y0[i % (n / 2)] + w[i] * y1[i % (n / 2)];
	}

	return y;
}

complexD MFCC_Class::getWi(int i, int n)
{
	double alpha = 2 * M_PI * i / n;
	return complexD(std::cos(alpha), std::sin(alpha));
}

std::vector<double> MFCC_Class::windowHamming(std::vector<complexD> frame)
{
	int n = frame.size();
	std::vector<double> frameD(n);
	double h;
	for (int i = 0; i < n; i++)
	{
		h = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));
		frame[i] *= h;
		frameD[i] = std::sqrt(std::norm(frame[i]));
	}
	return frameD;
}

std::vector<std::vector<double> > MFCC_Class::getMelFilters(int mfccLength, int frameLength, int frequency, int frequencyMin, int frequencyMax)
{
	std::vector<double> points(mfccLength + 2);
	points[0] = convertToMel(frequencyMin);
	points[mfccLength + 1] = frequencyMax;

	double step = (convertToMel(points[mfccLength + 1]) - points[0]) / (mfccLength + 1);
	for (int i = 1; i < mfccLength + 1; i++)
	{
		points[i] = points[i - 1] + step;
		points[i - 1] = superpose(convertFromMel(points[i - 1]), frameLength, frequency);
	}
	points[mfccLength] = superpose(convertFromMel(points[mfccLength]), frameLength, frequency);
	points[mfccLength + 1] = superpose(points[mfccLength + 1], frameLength, frequency);

	std::vector<std::vector<double> > filterBanks(mfccLength);
	for (int i = 0; i < mfccLength; i++)
	{
		filterBanks[i].resize(frameLength);
	}

	for (int i = 1; i < mfccLength + 1; i++)
	{
		for (int j = 0; j < frameLength; j++) {

			if (points[i - 1] <= j && j <= points[i])
			{
				filterBanks[i - 1][j] = (j - points[i - 1]) / (points[i] - points[i - 1]);
			}
			else if (points[i] < j && j <= points[i + 1])
			{
				filterBanks[i - 1][j] = (points[i + 1] - j) / (points[i + 1] - points[i]);
			}
			else
			{
				filterBanks[i - 1][j] = 0;
			}
		}
	}

	return filterBanks;
}

std::vector<double> MFCC_Class::logEnergySpectrum(std::vector<double> spectrum, std::vector<std::vector<double> > melFilters, int mfccCount) {

	std::vector<double> logPower(mfccCount);

	for (int i = 0; i < mfccCount; i++)
	{
		logPower[i] = 0.;

		for (int j = 0; j < spectrum.size(); j++)
		{
			logPower[i] += melFilters[i][j] * pow(spectrum[j], 2);
		}
		//??
		logPower[i] = std::log(logPower[i]);
	}

	return logPower;
}

std::vector<double> MFCC_Class::DCT(std::vector<double> data)
{
	int length = data.size();
	std::vector<double> cepstral(length);

	for (int i = 0; i < length; i++)
	{
		cepstral[i] = 0;
		for (int j = 0; j < length; j++)
		{
			cepstral[i] += data[j] * cos(M_PI * i * (j + 1.0 / 2.0) / length);
		}
	}

	return cepstral;
}

std::vector<double> MFCC_Class::computeMFCC(std::vector<double> frame, int frequency, int mfccCount, int frequencyMin, int frequencyMax)
{
	int frameLength = frame.size();
	std::vector<complexD> frameC(frameLength);

	for (int i = 0; i < frameLength; i++)
	{
		frameC[i] = complexD(frame[i], 0);
	}

	int n = (int)std::ceil(std::log2(frameLength));
	int maxLenght = (int)std::pow(2, n);
	frameC.resize(maxLenght, complexD(0, 0));

	frameC = FFT(frameC);
	//frameC.resize(frameLength);
	std::vector<double> spectrum = windowHamming(frameC);

	frameLength = maxLenght;
	std::vector<std::vector<double> > melFilters = getMelFilters(mfccCount, frameLength, frequency, frequencyMin, frequencyMax);
	std::vector<double> logPower = logEnergySpectrum(spectrum, melFilters, mfccCount);
	std::vector<double> mfcc = DCT(logPower);

	return mfcc;
}

