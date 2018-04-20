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
	if (m_verbose)
	{

	}
}


int abc(int argc, char* argv[])
{




	constexpr unsigned int Dimension = 3;
	constexpr unsigned int numOfBayesianClasses = 2;
	float variance = 0.3;
	float weight = 0.5;
	using ImageType = itk::Image<float, Dimension>;
	using LabelImageType = itk::Image<unsigned char, Dimension>;

	using ImageReaderType = itk::ImageFileReader<ImageType>;
	ImageReaderType::Pointer imageReader = ImageReaderType::New();
	//imageReader->SetFileName("D:/projects/itkBayesianFilter/data/BrainT1Slice.nii.gz");
	//imageReader->SetFileName("D:/projects/itkBayesianFilter/data/image.nii.gz");
	imageReader->SetFileName("D:/projects/itkBayesianFilter/data/teeth/image.nii");
	imageReader->Update();

	using LabelImageReader = itk::ImageFileReader<LabelImageType>;
	LabelImageReader::Pointer labelReader = LabelImageReader::New();
	//labelReader->SetFileName("D:/projects/itkBayesianFilter/data/BrainT1Slice_labeled.nii.gz");
	//labelReader->SetFileName("D:/projects/itkBayesianFilter/data/label_teeth_denoise.nii.gz");
	labelReader->SetFileName("D:/projects/itkBayesianFilter/data/teeth/label_teeth_IL.nii.gz");
	labelReader->Update();

	std::cout << "finish reading image..." << std::endl;

	using BayesianInitializerType = itk::BayesianClassifierInitializationImageFilter<ImageType>;
	BayesianInitializerType::Pointer bayesianInitializer = BayesianInitializerType::New();
	bayesianInitializer->SetInput(imageReader->GetOutput());
	bayesianInitializer->SetNumberOfClasses(numOfBayesianClasses);
	bayesianInitializer->Update();

	std::cout << "finish perfoming bayesian initialization..." << std::endl;

	// extract membership image
	using MembershipImageConstantIteratorType = itk::ImageRegionConstIterator< BayesianInitializerType::OutputImageType>;
	using MembershipImageIteratorType = itk::ImageRegionIterator<ImageType>;
	MembershipImageConstantIteratorType bayesianInitializeIterator(bayesianInitializer->GetOutput(),
		bayesianInitializer->GetOutput()->GetLargestPossibleRegion());

	std::vector<ImageType::Pointer> bayesianInitializerMembershipImageList;
	for (int i = 0; i < numOfBayesianClasses; i++)
	{
		// allocate the membership image
		ImageType::Pointer membershipImage = ImageType::New();
		bayesianInitializerMembershipImageList.push_back(membershipImage);
		membershipImage->CopyInformation(bayesianInitializer->GetOutput());
		membershipImage->SetBufferedRegion(bayesianInitializer->GetOutput()->GetBufferedRegion());
		membershipImage->SetRequestedRegion(bayesianInitializer->GetOutput()->GetRequestedRegion());
		membershipImage->Allocate();


		MembershipImageIteratorType membershipIterator(membershipImage, membershipImage->GetBufferedRegion());

		// iterate through the membership image
		bayesianInitializeIterator.GoToBegin();
		membershipIterator.GoToBegin();

		while (!bayesianInitializeIterator.IsAtEnd())
		{
			membershipIterator.Set(bayesianInitializeIterator.Get()[i]);
			++membershipIterator;
			++bayesianInitializeIterator;
		}
	}

	// bayesian classification
	using PriorType = float;
	using PosteriorType = float;
	using ClassifierFilterType = itk::BayesianClassifierImageFilter<
		BayesianInitializerType::OutputImageType, LabelImageType::PixelType,
		PosteriorType, PriorType >;
	ClassifierFilterType::Pointer bayesianClassifierFilter = ClassifierFilterType::New();
	bayesianClassifierFilter->SetInput(bayesianInitializer->GetOutput());

	////filter->SetNumberOfSmoothingIterations(1);
	////using ExtractedComponentImageType = ClassifierFilterType::ExtractedComponentImageType;
	////using SmoothingFilterType = itk::GradientAnisotropicDiffusionImageFilter<ExtractedComponentImageType, ExtractedComponentImageType >;
	////SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
	////smoother->SetNumberOfIterations(5);
	////smoother->SetTimeStep(0.125);
	////smoother->SetConductanceParameter(3);
	////filter->SetSmoothingFilter(smoother);
	//bayesianClassifierFilter->Update();

	//std::cout << "finish bayesian classification..." << std::endl;

	// calculate mean in each class with bayesian segmentation output
	using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter<ImageType, LabelImageType>;
	LabelStatisticsImageFilterType::Pointer bayesianOutputStatisticFilter = LabelStatisticsImageFilterType::New();
	bayesianOutputStatisticFilter->SetInput(imageReader->GetOutput());
	bayesianOutputStatisticFilter->SetLabelInput(bayesianClassifierFilter->GetOutput());
	bayesianOutputStatisticFilter->Update();

	// calculate the mean in each class with user input
	LabelStatisticsImageFilterType::Pointer userLabelStatisticFilter = LabelStatisticsImageFilterType::New();
	userLabelStatisticFilter->SetInput(imageReader->GetOutput());
	userLabelStatisticFilter->SetLabelInput(labelReader->GetOutput());
	userLabelStatisticFilter->Update();

	// create gaussian mask image
	using GaussianMembershipImageType = itk::VectorImage<float, Dimension>;
	GaussianMembershipImageType::Pointer gaussianMemebership = GaussianMembershipImageType::New();
	GaussianMembershipImageType::IndexType start;
	GaussianMembershipImageType::SizeType size;
	for (int i = 0; i < Dimension; i++)
	{
		start.SetElement(i, labelReader->GetOutput()->GetBufferedRegion().GetIndex()[i]);
		size.SetElement(i, labelReader->GetOutput()->GetBufferedRegion().GetSize()[i]);
	}
	GaussianMembershipImageType::RegionType region(start, size);
	gaussianMemebership->SetRegions(region);
	gaussianMemebership->SetNumberOfComponentsPerPixel(userLabelStatisticFilter->GetNumberOfLabels());
	gaussianMemebership->Allocate();

	// gaussian weighted bayesian initialization membership 
	int count = 0;

	for (auto vIt = userLabelStatisticFilter->GetValidLabelValues().begin();
		vIt != userLabelStatisticFilter->GetValidLabelValues().end();
		++vIt)
	{
		if (userLabelStatisticFilter->HasLabel(*vIt))
		{
			unsigned short labelValue = *vIt;
			std::cout << "label: " << labelValue << std::endl;
			std::cout << "mean: " << userLabelStatisticFilter->GetMean(labelValue) << std::endl;
			std::cout << "variance: " << userLabelStatisticFilter->GetVariance(labelValue) << std::endl;

			// find the bayesian segmented class with nearest mean
			float diff = std::numeric_limits<float>::max();
			unsigned short correspondBayesianClass = 0;
			for (auto bayesianLabelStatIt = bayesianOutputStatisticFilter->GetValidLabelValues().begin();
				bayesianLabelStatIt != bayesianOutputStatisticFilter->GetValidLabelValues().end();
				++bayesianLabelStatIt)
			{
				if (bayesianOutputStatisticFilter->HasLabel(*bayesianLabelStatIt))
				{
					if (std::abs(bayesianOutputStatisticFilter->GetMean(*bayesianLabelStatIt) - userLabelStatisticFilter->GetMean(labelValue)) < diff)
					{
						diff = std::abs(bayesianOutputStatisticFilter->GetMean(*bayesianLabelStatIt) - userLabelStatisticFilter->GetMean(labelValue));
						correspondBayesianClass = *bayesianLabelStatIt;
					}
				}
			}

			std::cout << "corresponding bayesian class: " << correspondBayesianClass << std::endl;
			std::cout << "mean: " << bayesianOutputStatisticFilter->GetMean(correspondBayesianClass) << std::endl;

			// calculate statistic of each membership class
			using StatisticsImageFilterType = itk::StatisticsImageFilter<ImageType>;
			StatisticsImageFilterType::Pointer bayesianMembershipStatFilter = StatisticsImageFilterType::New();
			bayesianMembershipStatFilter->SetInput(bayesianInitializerMembershipImageList.at(correspondBayesianClass));
			bayesianMembershipStatFilter->Update();

			std::cout << "Bayesian membership class: " << correspondBayesianClass << std::endl;
			std::cout << "Mean: " << bayesianMembershipStatFilter->GetMean() << std::endl;
			std::cout << "Std.: " << bayesianMembershipStatFilter->GetSigma() << std::endl;


			std::cout << "***************************************" << std::endl;

			// extract the label
			using BinaryThresholdImageFilterType = itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>;
			BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
			thresholdFilter->SetInput(labelReader->GetOutput());
			thresholdFilter->SetUpperThreshold(labelValue);
			thresholdFilter->SetLowerThreshold(labelValue);
			thresholdFilter->SetInsideValue(1);
			thresholdFilter->SetOutsideValue(0);
			thresholdFilter->Update();

			// apply gaussian blur on the label
			using GaussianFilterType = itk::DiscreteGaussianImageFilter<LabelImageType, ImageType>;
			GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
			gaussianFilter->SetInput(thresholdFilter->GetOutput());
			gaussianFilter->SetVariance(variance);
			gaussianFilter->Update();

			using ConstLabelImageIteratorType = itk::ImageRegionConstIterator< ImageType >;
			using ConstMembershipIteratorType = itk::ImageRegionConstIterator<BayesianInitializerType::OutputImageType>;
			using GaussianMembershipIteratorType = itk::ImageRegionIterator<GaussianMembershipImageType>;
			ConstLabelImageIteratorType labelIterator(gaussianFilter->GetOutput(),
				gaussianFilter->GetOutput()->GetLargestPossibleRegion());
			ConstMembershipIteratorType bayesianInitializerIterator(bayesianInitializer->GetOutput(),
				bayesianInitializer->GetOutput()->GetLargestPossibleRegion());
			GaussianMembershipIteratorType gaussianMembershipIterator(gaussianMemebership,
				gaussianMemebership->GetLargestPossibleRegion());

			labelIterator.GoToBegin();
			bayesianInitializerIterator.GoToBegin();
			gaussianMembershipIterator.GoToBegin();

			while (!labelIterator.IsAtEnd())
			{
				gaussianMembershipIterator.Get().SetElement(count, ((1.0 - weight)*bayesianInitializerIterator.Get()[correspondBayesianClass] * 0.5 + weight*labelIterator.Get() * bayesianMembershipStatFilter->GetMean()) / ((0.5 + bayesianMembershipStatFilter->GetMean()) / 2.0));
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

			count++;
		}
	}

	//// save the final membership image
	//for (int i = 0; i < gaussianMemebership->GetNumberOfComponentsPerPixel(); i++)
	//{
	//	using ExtractedComponentImageType = itk::Image< float, Dimension >;

	//	ExtractedComponentImageType::Pointer extractedComponentImage = ExtractedComponentImageType::New();
	//	extractedComponentImage->CopyInformation(gaussianMemebership);
	//	extractedComponentImage->SetBufferedRegion(gaussianMemebership->GetBufferedRegion());
	//	extractedComponentImage->SetRequestedRegion(gaussianMemebership->GetRequestedRegion());
	//	extractedComponentImage->Allocate();
	//	using ConstIteratorType = itk::ImageRegionConstIterator< GaussianMembershipImageType >;
	//	using IteratorType = itk::ImageRegionIterator< ExtractedComponentImageType >;
	//	ConstIteratorType cit(gaussianMemebership,
	//		gaussianMemebership->GetBufferedRegion());
	//	IteratorType it(extractedComponentImage,
	//		extractedComponentImage->GetLargestPossibleRegion());

	//	cit.GoToBegin();
	//	it.GoToBegin();

	//	while (!cit.IsAtEnd())
	//	{
	//		it.Set(cit.Get()[i]);

	//		++it;
	//		++cit;
	//	}


	//	using WriterType = itk::ImageFileWriter< ExtractedComponentImageType >;
	//	WriterType::Pointer membershipWriter = WriterType::New();
	//	membershipWriter->SetFileName("D:/projects/itkBayesianFilter/data/membership_" + std::to_string(i)+ ".nii.gz");
	//	membershipWriter->SetInput(extractedComponentImage);
	//	membershipWriter->Update();
	//}

	// final bayesian classification
	ClassifierFilterType::Pointer bayesianClassifierFilter2 = ClassifierFilterType::New();
	bayesianClassifierFilter2->SetInput(gaussianMemebership);
	bayesianClassifierFilter2->Update();

	//using OutputImageType = itk::Image< unsigned char, Dimension >;
	using LabelWriterType = itk::ImageFileWriter<LabelImageType>;
	LabelWriterType::Pointer labelWriter = LabelWriterType::New();
	labelWriter->SetInput(bayesianClassifierFilter2->GetOutput());
	//writer->SetFileName("D:/projects/itkBayesianFilter/data/class_" + std::to_string(i) + ".nii");
	//labelWriter->SetFileName("D:/projects/itkBayesianFilter/data/label_bayesian_new.nii");
	labelWriter->SetFileName("D:/projects/itkBayesianFilter/data/teeth/label_teeth_IL_bayesian.nii.gz");

	labelWriter->Update();
}

