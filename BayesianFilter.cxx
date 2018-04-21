#include "BayesianFilter.h"


#include "itkImage.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkStatisticsImageFilter.h"

BayesianFilter::BayesianFilter()
{
	// parameter initialization
	m_numOfInitClasses = 2;
	m_variance = 0.3;
	m_weight = 0.5;
	m_verbose = true;
	m_smooth = false;
	m_smoothItr = 1;
	m_gadItr = 5;
	m_timeStep = 0.125;
	m_conductance = 3;
	m_output = LabelImageType::New();
}

BayesianFilter::~BayesianFilter()
{
}

void BayesianFilter::SetNumberOfBayesianInitialClasses(unsigned int numOfInitClasses)
{
	m_numOfInitClasses = numOfInitClasses;
}

void BayesianFilter::SetGaussianBlurVariance(float variance)
{
	m_variance = variance;
}

void BayesianFilter::SetLabelWeight(float weight)
{
	m_weight = weight;
}

void BayesianFilter::SetVerbose(bool verbose)
{
	m_verbose = verbose;
}

void BayesianFilter::SetSmooth(bool smooth)
{
	m_smooth = smooth;
}

void BayesianFilter::SetNumberOfSmoothingIterations(unsigned int smoothItr)
{
	m_smoothItr = smoothItr;
}

void BayesianFilter::SetNumberOfGradientAnistropicDiffusionIterations(unsigned int gadItr)
{
	m_gadItr = gadItr;
}

void BayesianFilter::SetSmoothingTimeStep(float timeStep)
{
	m_timeStep = timeStep;
}

void BayesianFilter::SetConductanceParameter(float conductance)
{
	m_conductance = conductance;
}

void BayesianFilter::SetImage(ImageType::Pointer inputImage)
{
	m_inputImage->Graft(inputImage);
}

void BayesianFilter::SetLabel(LabelImageType::Pointer labelImage)
{
	m_labelImage->Graft(labelImage);
}

void BayesianFilter::Run()
{
	// check input is valid
	if (m_inputImage == nullptr || m_labelImage == nullptr)
	{
		std::cerr << "Invalid input data" << std::endl;
		return;
	}

	// check input image and label have same size
	if (m_inputImage->GetBufferedRegion().GetSize()[0] != m_labelImage->GetBufferedRegion().GetSize()[0] ||
		m_inputImage->GetBufferedRegion().GetSize()[1] != m_labelImage->GetBufferedRegion().GetSize()[1] || 
		m_inputImage->GetBufferedRegion().GetSize()[2] != m_labelImage->GetBufferedRegion().GetSize()[2])
	{
		std::cerr << "Input image and label should have same size" << std::endl;
		return;
	}

	// Bayesian initialization
	if (m_verbose)
	{
		std::cout << "Bayesian refinement in progress..." << std::endl;
		std::cout << "Start Bayesian initialization..." << std::endl;
	}

	BayesianInitializerType::Pointer bayesianInitializer = BayesianInitializerType::New();
	bayesianInitializer->SetInput(m_inputImage);
	bayesianInitializer->SetNumberOfClasses(m_numOfInitClasses);
	bayesianInitializer->Update();
	
	// Extract membership image
	if (m_verbose)
	{
		std::cout << "Extracting membership images from Bayesian intialization..." << std::endl;
	}

	MembershipImageRegionConstIteratorType bayesianInitializerIterator(bayesianInitializer->GetOutput(),
		bayesianInitializer->GetOutput()->GetBufferedRegion());
	std::vector<ImageType::Pointer> bayesianInitializerMembershipImageList;
	for (int i = 0;i < m_numOfInitClasses; i++)
	{
		// allocate the membership image of each initialized class
		ImageType::Pointer membershipImage = ImageType::New();
		bayesianInitializerMembershipImageList.push_back(membershipImage);
		membershipImage->CopyInformation(bayesianInitializer->GetOutput());
		membershipImage->SetBufferedRegion(bayesianInitializer->GetOutput()->GetBufferedRegion());
		membershipImage->SetRequestedRegion(bayesianInitializer->GetOutput()->GetRequestedRegion());
		membershipImage->Allocate();

		// copy memberhsip image from Bayesian initializer output
		ImageRegionIteratorType membershipIterator(membershipImage,
			membershipImage->GetBufferedRegion());
		bayesianInitializerIterator.GoToBegin();
		membershipIterator.GoToBegin();

		while (!bayesianInitializerIterator.IsAtEnd())
		{
			membershipIterator.Set(bayesianInitializerIterator.Get()[i]);
			++membershipIterator;
			++bayesianInitializerIterator;
		}
	}

	// Bayesian classification using Kmean
	if (m_verbose)
	{
		std::cout << "Performing Bayesian classification with Kmean membership..." << std::endl;
	}

	// bayesian classification
	BayesianClassifierFilterType::Pointer bayesianClassifierFiler = BayesianClassifierFilterType::New();
	bayesianClassifierFiler->SetInput(bayesianInitializer->GetOutput());

	// apply smoothing on the Bayesian classification
	if (m_smooth)
	{
		bayesianClassifierFiler->SetNumberOfSmoothingIterations(1);
		GADSmoothingFilterType::Pointer smoother = GADSmoothingFilterType::New();
		smoother->SetNumberOfIterations(5);
		smoother->SetTimeStep(0.125);
		smoother->SetConductanceParameter(3);
		bayesianClassifierFiler->SetSmoothingFilter(smoother);
	}
	bayesianClassifierFiler->Update();

	// calculate mean in each class with bayesian segmentation output
	if (m_verbose)
	{
		std::cout << "Calculating statistics on segmentation labels..." << std::endl;
	}

	LabelStatisticsImageFilterType::Pointer bayesianLabelStatisticFilter = LabelStatisticsImageFilterType::New();
	bayesianLabelStatisticFilter->SetInput(m_inputImage);
	bayesianLabelStatisticFilter->SetLabelInput(bayesianClassifierFiler->GetOutput());
	bayesianLabelStatisticFilter->Update();

	// calculate the mean in each class with user input

	LabelStatisticsImageFilterType::Pointer userLabelStatisticFilter = LabelStatisticsImageFilterType::New();
	userLabelStatisticFilter->SetInput(m_inputImage);
	userLabelStatisticFilter->SetLabelInput(m_labelImage);
	userLabelStatisticFilter->Update();

	// gaussian weighted bayesian initialization membership 
	if (m_verbose)
	{
		std::cout << "Calculating user label weighted membership..." << std::endl;
	}

	// create gaussian mask image
	StatisticsLabelImageFilterType::Pointer labelStatFilter = StatisticsLabelImageFilterType::New();
	labelStatFilter->SetInput(m_labelImage);
	labelStatFilter->Update();

	MembershipImageType::Pointer gaussianMemebership = MembershipImageType::New();
	MembershipImageType::IndexType start;
	MembershipImageType::SizeType size;
	for (int i = 0; i < 3; i++)
	{
		start.SetElement(i, m_labelImage->GetBufferedRegion().GetIndex()[i]);
		size.SetElement(i, m_labelImage->GetBufferedRegion().GetSize()[i]);
	}
	MembershipImageType::RegionType region(start, size);
	gaussianMemebership->SetRegions(region);
	gaussianMemebership->SetNumberOfComponentsPerPixel(labelStatFilter->GetMaximum()+1); // note that there is 0 class
	gaussianMemebership->Allocate();

	// calculate the membership for each class in user label input
	for (auto vIt = userLabelStatisticFilter->GetValidLabelValues().begin();
		vIt != userLabelStatisticFilter->GetValidLabelValues().end();
		++vIt)
	{
		if (userLabelStatisticFilter->HasLabel(*vIt))
		{
			unsigned short labelValue = *vIt;
			//std::cout << "label: " << labelValue << std::endl;
			//std::cout << "mean: " << userLabelStatisticFilter->GetMean(labelValue) << std::endl;

			// find the bayesian segmented class with nearest mean
			float diff = std::numeric_limits<float>::max();
			unsigned short correspondBayesianClass = 0;
			for (auto bayesianLabelStatIt = bayesianLabelStatisticFilter->GetValidLabelValues().begin();
				bayesianLabelStatIt != bayesianLabelStatisticFilter->GetValidLabelValues().end();
				++bayesianLabelStatIt)
			{
				if (bayesianLabelStatisticFilter->HasLabel(*bayesianLabelStatIt))
				{
					if (std::abs(bayesianLabelStatisticFilter->GetMean(*bayesianLabelStatIt) - userLabelStatisticFilter->GetMean(labelValue)) < diff)
					{
						diff = std::abs(bayesianLabelStatisticFilter->GetMean(*bayesianLabelStatIt) - userLabelStatisticFilter->GetMean(labelValue));
						correspondBayesianClass = *bayesianLabelStatIt;
					}
				}
			}

			//std::cout << "corresponding bayesian class: " << correspondBayesianClass << std::endl;
			//std::cout << "mean: " << bayesianLabelStatisticFilter->GetMean(correspondBayesianClass) << std::endl;

			// calculate statistic of each membership class
			StatisticsImageFilterType::Pointer bayesianMembershipStatFilter = StatisticsImageFilterType::New();
			bayesianMembershipStatFilter->SetInput(bayesianInitializerMembershipImageList.at(correspondBayesianClass));
			bayesianMembershipStatFilter->Update();

			//std::cout << "Bayesian membership class: " << correspondBayesianClass << std::endl;
			//std::cout << "Mean: " << bayesianMembershipStatFilter->GetMean() << std::endl;

			//std::cout << "***************************************" << std::endl;

			// extract the label
			using BinaryThresholdImageFilterType = itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>;
			BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
			thresholdFilter->SetInput(m_labelImage);
			thresholdFilter->SetUpperThreshold(labelValue);
			thresholdFilter->SetLowerThreshold(labelValue);
			thresholdFilter->SetInsideValue(1);
			thresholdFilter->SetOutsideValue(0);
			thresholdFilter->Update();

			// apply gaussian blur on the label
			GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
			gaussianFilter->SetInput(thresholdFilter->GetOutput());
			gaussianFilter->SetVariance(m_variance);
			gaussianFilter->Update();

			ImageRegionConstIteratorType labelIterator(gaussianFilter->GetOutput(),
				gaussianFilter->GetOutput()->GetLargestPossibleRegion());
			//MembershipImageRegionConstIteratorType bayesianInitializerIterator(bayesianInitializer->GetOutput(),
			//	bayesianInitializer->GetOutput()->GetLargestPossibleRegion());
			MembershipImageRegionIteratorType gaussianMembershipIterator(gaussianMemebership,
				gaussianMemebership->GetLargestPossibleRegion());

			labelIterator.GoToBegin();
			bayesianInitializerIterator.GoToBegin();
			gaussianMembershipIterator.GoToBegin();

			while (!labelIterator.IsAtEnd())
			{
				gaussianMembershipIterator.Get().SetElement(labelValue, ((1.0 - m_weight)*bayesianInitializerIterator.Get()[correspondBayesianClass] * 0.5 + m_weight*labelIterator.Get() * bayesianMembershipStatFilter->GetMean()) / ((0.5 + bayesianMembershipStatFilter->GetMean()) / 2.0));
				//if (gaussianMembershipIterator.Get().GetElement(count) < 0)
				//{
				//	std::cout << "correspondBayesianClass: " << correspondBayesianClass << std::endl;
				//	std::cout << "bayesianInitializeIterator.Get()[correspondBayesianClass]: " << bayesianInitializerIterator.Get()[correspondBayesianClass] << std::endl;
				//	std::cout << "labelIterator.Get(): " << labelIterator.Get() << std::endl;
				//	std::cout << "bayesianMembershipStatFilter->GetMean(): " << bayesianMembershipStatFilter->GetMean() << std::endl;
				//}
				++gaussianMembershipIterator;
				++bayesianInitializerIterator;
				++labelIterator;
			}
		}
	}

	// final Bayesian segmentation
	if (m_verbose)
	{
		std::cout << "Performing label refinement..." << std::endl;
	}
	// final bayesian classification
	BayesianClassifierFilterType::Pointer bayesianClassifierFilter2 = BayesianClassifierFilterType::New();
	bayesianClassifierFilter2->SetInput(gaussianMemebership);
	bayesianClassifierFilter2->Update();

	m_output->Graft(bayesianClassifierFilter2->GetOutput());
}

LabelImageType::Pointer BayesianFilter::GetOutput()
{
	return m_output;
}

void BayesianFilter::SetOutput(LabelImageType::Pointer output)
{
	m_output = output;
}
