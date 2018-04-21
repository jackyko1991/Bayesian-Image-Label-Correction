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
	void SetSmooth(bool smooth);
	void SetNumberOfSmoothingIterations(unsigned int smoothItr);
	void SetNumberOfGradientAnistropicDiffusionIterations(unsigned int gadItr);
	void SetSmoothingTimeStep(float timeStep);
	void SetConductanceParameter(float conductance);
	void Run();
	LabelImageType::Pointer GetOutput();
	void SetOutput(LabelImageType::Pointer output);

private:
	ImageType::Pointer m_inputImage;
	LabelImageType::Pointer m_labelImage;
	LabelImageType::Pointer m_output;
	unsigned int m_numOfInitClasses;
	float m_variance;
	float m_weight;
	bool m_verbose;
	bool m_smooth;
	unsigned int m_smoothItr;
	unsigned int m_gadItr;
	float m_timeStep;
	float m_conductance;
};

#endif