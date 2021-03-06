///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            *** Model.hpp ***                              //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Model_hpp
#define Model_hpp

#include <vector>
#include <phynet/Network.hpp>
#include <phynet/Loss.hpp>
#include <phynet/Optimizer.hpp>
#include <gendat/Operators.hpp>
#include <chrono>

template <typename T>
class Model  
{
	public:
		Model(const std::vector<Network<T>>& networks, const Loss<T>& loss, 
			  const Optimizer<T>& optimizer);

		void learn_from(const Dataset<T>& dataset); 
		void train(const Dataset<T>& dataset);
		void reset(void);

		T mse(const Dataset<T>& dataset);

		void append_metrics(const Dataset<T> &dataset, int trial, int epoch, std::string fpath);

		void write_overlap(const Dataset<T> &dataset, std::string fpath);
		void write_magnetization(const Dataset<T> &dataset, std::string fpath);

		void write_radial_visualization(const Dataset<T> &dataset, std::string fpath);
		void write_entanglement_entropy(const Dataset<T> &dataset, std::string fpath);

		void print_average_overlap(const Dataset<T>& dataset);
		void print_average_sz_error(const Dataset<T> &dataset);
		void print_inference_time(const Dataset<T> &dataset);

		void append_wandb_for_radviz(std::string fpath);

		//void read(std::string fpath);
		//void save(std::string fpath);

		//void write_coefficients(const Dataset<T> &dataset, std::string fpath, int epoch);
		//void write_schrodinger_error(const Dataset<T> &dataset, std::string fpath, int epoch);
		void write_lyapunov_estimate(const Dataset<T> &dataset, std::string fpath);

		T lyapunov_estimate(const Dataset<T> &dataset);
		T overlap(const Dataset<T> &dataset);

	private:
		T entanglement_entropy(Eigen::Matrix<T, Eigen::Dynamic, 1> psi) const;

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			reduced_density_matrix(Eigen::Matrix<T, Eigen::Dynamic, 1> psi) const;

		std::vector<Network<T>> networks;
		Loss<T> loss; 
		Optimizer<T> optimizer;		

		void pretty_print(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
};
	
#endif /* Model_hpp */

