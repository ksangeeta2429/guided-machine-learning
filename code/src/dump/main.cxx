///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** dump main.cxx ***                           //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <dump/Reader.hpp>

int main(int argc, char *argv[])
{
	if (argc != 3 && argc != 4)
	{
		if (argv[0])
		{
			std::cout << " Usage: " << argv[0];
			std::cout << " <num_qubits> ";
			std::cout << " <fpath>";
			std::cout << " [-csv]";
			std::cout << std::endl;
		}
		else
		{
			std::cout << " Usage: dump";
			std::cout << " <num_qubits> ";
			std::cout << " <fpath>";
			std::cout << " [-csv]";
			std::cout << std::endl;
		}

		exit(-1);
	}

	std::stringstream a1(argv[1]);

	int num_qubits;

	a1 >> num_qubits;

	Reader<double> reader(num_qubits, argv[2]);	

	reader.read();

	reader.print();

	if (argc == 4) 
	{
		reader.write_csv();
		std::cout << "Fields written to input.csv, wavefx's to output.csv\n";
	}

	return 0;
}


