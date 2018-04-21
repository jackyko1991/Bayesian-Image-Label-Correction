#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkBayesianClassifierInitializationImageFilter.h>
#include <itkBayesianClassifierImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>

typedef itk::Image<float, 3> ImageType;
typedef itk::Image<unsigned short, 3> LabelImageType;
typedef itk::VectorImage<float, 3> MembershipImageType;

typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

typedef itk::ImageFileWriter<MembershipImageType> MemebershipImageWriterType;
typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;

typedef itk::BayesianClassifierInitializationImageFilter<ImageType> BayesianInitializerType;
typedef itk::BayesianClassifierImageFilter<MembershipImageType, LabelImageType::PixelType, float, float> BayesianClassifierFilterType;

typedef itk::ImageRegionConstIterator<MembershipImageType> MembershipImageRegionConstIteratorType;
typedef itk::ImageRegionIterator<ImageType> ImageRegionIteratorType;
typedef itk::ImageRegionIterator<MembershipImageType> MembershipImageRegionIteratorType;
typedef itk::ImageRegionConstIterator<ImageType> ImageRegionConstIteratorType;

typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> GADSmoothingFilterType;
typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> LabelStatisticsImageFilterType;
typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
typedef itk::StatisticsImageFilter<LabelImageType> StatisticsLabelImageFilterType;
typedef itk::DiscreteGaussianImageFilter<LabelImageType, ImageType> GaussianFilterType;

#endif