#pragma once

// class to keep parameters neatly organized



// C++ STL
#include <vector>								// std::vector
#include <random>								// default_random_engine, normal_distribution<double>, uniform_real_distribution


/// @brief replace integers with verbose names
enum class model : unsigned int {
	spherical = 0, 
	confined_3d = 1, 
	confined_3d_torqueNoise = 2, 
	bistable_potential = 3,
	confined_3d_torqueNoise_telegraph = 4,
	planar_2d_torqueNoise_telegraph = 5,
	planar_2d_ornstein_uhlenbeck = 6,
	planar_1d_ornstein_uhlenbeck = 7,
	spherical_gravity_liquid_crystal = 8
};

/// @brief simulation parameters used throughout the program
class Params {
	
	public:
		
		// methods
		Params(const int argc, char *argv[]);
		
		// attributes
		double diffRot, diffS;							// 
		double rDroplet, rForce;						// droplet radius, droplet radius + bacterium radius
		double dragFudgeFactor;							// drag Fudge Factor
		double swimFactor, gravFactor, lcFactor;
		double dt;										// time step
		double timeFactor, c;
		double rateRT, rateTR;							// transition rates in telegraph process
		double k;										// spring constant in Ornstein Uhlenbeck process
		int nStates, nSaveStates;						// 
		int nSystems, nVariables, nHistoryStates;		// 
		model geometryChoice;							// which model to use (int as named enum class)
		int noiseSeed;									// noise seed
		bool hetQ;										// use heterogeneity?
		std::vector<double> het;
		std::normal_distribution<double> normal_distro;
		std::uniform_real_distribution<double> uniform_distro;
		std::mt19937_64 generator;								// random number generator
		std::string pthout;
		
	private:
		
		// methods
		void readCmdline(int argc, char *argv[]);
		void create_heterogeneity_values(const double scalingFactor);
		void write_heterogeneity_values();
		
		// attributes
		std::vector<double> cmdParams;
};
