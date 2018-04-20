#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

typedef itk::Image<float, 3> ImageType;
typedef itk::Image<unsigned short, 3> LabelImageType;
typedef itk::VectorImage<float, 3> MembershipImageType;


typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

typedef itk::ImageFileWriter<MembershipImageType> MemebershipImageWriterType;
typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;

#endif