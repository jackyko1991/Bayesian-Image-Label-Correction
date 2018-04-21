#include "BayesianFilter.h"

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "iostream"

#include "typedef.h"
#include "itkImage.h"

int main(int argc, char* argv[])
{
	try
	{
		boost::program_options::options_description desc{
			"Adjust prelabeled images with Bayesian classification method." };
		desc.add_options()
			("help,h", "Print help")
			("input,i", boost::program_options::value<boost::filesystem::path>()->default_value("./data/image.png"), "Input image data")
			("label,l", boost::program_options::value<boost::filesystem::path>()->default_value("./data/label.png"), "Input label data")
			("output,o", boost::program_options::value<boost::filesystem::path>()->default_value("./data/output.png"), "Output label data")
			("classes,c", boost::program_options::value<unsigned int>()->default_value(2), "Number of classes used for Bayesian initialization")
			("variance,v", boost::program_options::value<float>()->default_value(0.3), "Variance of Gaussian blur")
			("weight,w", boost::program_options::value<float>()->default_value(0.5), "Weight of input label in membership")
			;

		boost::program_options::variables_map variableMap;
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), variableMap);

		if (variableMap.count("help") || argc == 1)
			std::cout << desc << std::endl;
		else
		{
			BayesianFilter bayesianFilter;

			// read input image
			if (variableMap.count("input"))
			{
				if (boost::filesystem::exists(variableMap["input"].as<boost::filesystem::path>()))
				{
					// read the image
					ImageReaderType::Pointer imageReader = ImageReaderType::New();
					imageReader->SetFileName(variableMap["input"].as<boost::filesystem::path>().string());

					try
					{
						imageReader->Update();
						bayesianFilter.SetImage(imageReader->GetOutput());
					}
					catch (itk::ExceptionObject & excp)
					{
						std::cerr << "Exception thrown " << std::endl;
						std::cerr << excp << std::endl;
						return EXIT_FAILURE;
					}
				}
				else
				{
					std::cerr << "Input image not exist!" << std::endl;
					return EXIT_FAILURE;
				}
			}

			// read input label
			if (variableMap.count("label"))
			{
				if (boost::filesystem::exists(variableMap["label"].as<boost::filesystem::path>()))
				{
					// read the image
					LabelReaderType::Pointer labelReader = LabelReaderType::New();
					labelReader->SetFileName(variableMap["label"].as<boost::filesystem::path>().string());

					try
					{
						labelReader->Update();
						bayesianFilter.SetLabel(labelReader->GetOutput());
					}
					catch (itk::ExceptionObject & excp)
					{
						std::cerr << "Exception thrown " << std::endl;
						std::cerr << excp << std::endl;
						return EXIT_FAILURE;
					}
				}
				else
				{
					std::cerr << "Input label not exist!" << std::endl;
					return EXIT_FAILURE;
				}
			}

			// set number of Bayesian initialization classes
			bayesianFilter.SetNumberOfBayesianInitialClasses(variableMap["classes"].as<unsigned int>());

			// set Gaussian blur variance
			if (variableMap["variance"].as<float>() >= 0)
			{
				bayesianFilter.SetNumberOfBayesianInitialClasses(variableMap["classes"].as<unsigned int>());
			}
			else
			{
				std::cerr << "Variance for Gaussian blur should be largeer than 0" << std::endl;
				return EXIT_FAILURE;
			}

			// run the code
			bayesianFilter.Run();

			// write the output
			LabelWriterType::Pointer labelWriter = LabelWriterType::New();
			labelWriter->SetFileName(variableMap["output"].as<boost::filesystem::path>().string());
			labelWriter->SetInput(bayesianFilter.GetOutput());
			labelWriter->Write();
		}
	}
	catch (const boost::program_options::error &ex)
	{
		std::cerr << ex.what() << std::endl;
		return 0;
	}
}