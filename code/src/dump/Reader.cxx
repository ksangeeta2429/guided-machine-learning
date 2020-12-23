///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Reader.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <dump/Reader.hpp>

template <typename T>
void Reader<T>::read(void)
{
	int dim = (int)(pow(2, num_qubits) + 0.5);	

	long offset = (3*num_qubits + dim + dim*dim)*sizeof(T);
	int num_instances = boost::filesystem::file_size(fpath)/offset;

	std::ifstream file(fpath.c_str(), std::ios::binary);

	if (!file.is_open())
	{
		std::cerr << "ERROR: COULD NOT OPEN " << fpath << '\n';
		exit(-1);
	}

	Fields<T> tmp_fields(num_qubits);	
	Eigen::Matrix<T, Eigen::Dynamic, 1> tmp_values(dim);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp_wavefx(dim, dim);
	
	if (file.eof())
	{
		std::cerr << "ERROR: FILESTREAM ERROR\n\n";
		exit(-1);
	}

	for (int i = 0; i < num_instances; ++i)
	{
		file.read((char*)tmp_fields.coupling.data(), num_qubits * sizeof(T));
		file.read((char*)tmp_fields.transverse.data(), num_qubits * sizeof(T));
		file.read((char*)tmp_fields.longitudinal.data(), num_qubits * sizeof(T));
		file.read((char*)tmp_values.data(), dim * sizeof(T));
		file.read((char*)tmp_wavefx.data(), dim * dim * sizeof(T));

		fields.push_back(tmp_fields);
		values.push_back(tmp_values);
		wavefx.push_back(tmp_wavefx);
	}
}

template <typename T>
void Reader<T>::write_csv(void) const 
{
	std::ofstream inputs("input.csv");
	std::ofstream outputs("output.csv");

	inputs << "J[1], Bx[1], Bz[1]... J[n], Bx[n], Bz[n]\n";	
	for (std::size_t i = 0; i < fields.size(); ++i)
	{
		inputs << std::setw(4) << i+1 << ',';
		for (std::size_t j = 0; j < fields[i].coupling.size(); ++j)
		{
			inputs << std::setw(10) << fields[i].coupling[j] << ',';
			inputs << std::setw(10) << fields[i].transverse[j] << ',';
			inputs << std::setw(10) << fields[i].longitudinal[j];

			if (j != fields[i].coupling.size() - 1) inputs << ',';
		}

		inputs << '\n';
	}

	outputs << "c[1], c[2], c[3]... c[n] \n";
	for (std::size_t i = 0; i < fields.size(); ++i)
	{
		outputs << std::setw(4) << i+1 << ',';
		for (int j = 0; j < wavefx[i].size(); ++j)
		{
			outputs << *(wavefx[i].data() + j);
			if (j != wavefx[i].size() - 1) outputs << ',';
		}
		outputs << '\n';
	}
}

template <typename T>
void Reader<T>::print(void) const
{
	for (std::size_t i = 0; i < fields.size(); ++i)
	{
		std::cout << std::scientific;
		std::cout << "--------------------------- INSTANCE " << i+1 << " ---------------------------\n";
		fields[i].print();

		std::cout << "Eigenvalues:\n";
		std::cout << values[i].transpose() << "\n\n";
		std::cout << "Eigenvectors:\n";
		std::cout << wavefx[i] << "\n";
		std::cout << "--------------------------------------------------------------------\n\n\n";
	}
}

template class Reader<double>;
