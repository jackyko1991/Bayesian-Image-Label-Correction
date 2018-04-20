#ifndef BAYESIANFILTER_H
#define BAYESIANFILTER_H

#include "typedef.h"
#include "itkImage.h"

class BayesianFilter
{
public:
	BayesianFilter();
	~BayesianFilter();
	void SetImage(ImageType::Pointer inputImage);
	void SetLabel(LabelImageType::Pointer labelImage);
	void SetNumberOfBayesianInitialClasses(unsigned int numOfInitClasses);
	void SetGaussianBlurVariance(float variance);
	void SetLabelWeight(float weight);
	void SetVerbose(bool verbose);
	void Run();

private:
	ImageType::Pointer m_inputImage;
	LabelImageType::Pointer m_labelImage;
	unsigned int m_numOfInitClasses;
	float m_variance;
	float m_weight;
	bool m_verbose;
};

#endif