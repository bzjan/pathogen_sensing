
// C++ STL
#include <iostream>								// std::cout, std::cerr, std::endl
#include <fstream>								// std::ofstream



// custom
#include "./Params.h"							// Params
#include "./utility_functions.h"				// get_executable_path(), print_std_vector(), blue,red,reset
#include "./droplet.h"							// calculate droplet parameters


/// @brief read parameters passed on commandline
/// @param argc number of arguments
/// @param argv values of arguments
void Params::readCmdline(const int argc, char *argv[]){
	
	// TODO: better cmdline parameter reader with named options
	const int nCmdParams = 14;
	this->cmdParams.resize(nCmdParams);
	
	if(argc == nCmdParams+1){					// read parameters from commandline, if correct number (+1 for name of program)
		for(int i=0; i<nCmdParams; i++) this->cmdParams[i] = std::stod(argv[i+1]);
	}else if(argc == 1){						// set default parameters if none are set
		this->cmdParams = {
			20.0e-6,					// e. coli speed u0Dim in m
			3.5,						// diffRot
			120.0,						// diffS
			0.2,						// timeFactor
			5.0e-6,						// rDropletDim in m
			0.0,						// asymmetry parameter in double well potential
			1.0,						// swimFactor
			100.0,						// tfinal in s
			1,							// number of systems solved in parallel
			0,							// noise seed
			0, 							// use heterogeneity; 0-false; 1-true
			1.1111, 					// transition rate from R -> T in 1/s; Lopez et al. 2019
			10.0, 						// transition rate from T -> R in 1/s; Lopez et al. 2019
			0.53 						// dragFudgeFactor
		};
	}else{								// user error
		printf("Error: Incorrect number (%d) of parameters are supplied!\n",argc); 
		for(int i=0; i<argc; i++) printf("%d) %s\n",i,argv[i]);
		exit(EXIT_FAILURE);
	}
	print_std_vector(this->cmdParams,"cmdParams");			// DEBUG
}


void Params::create_heterogeneity_values(const double scalingFactor){
	const double mean = this->swimFactor;
	// const double relativeWidth = 0.2;
	// const double stddev = relativeWidth*this->swimFactor;
	// const double u_stddev = 4.0e-6;  // speed is in (m /s)
	const double u_stddev = 3.0e-6;  // speed is in (m /s)
	const double stddev = u_stddev*scalingFactor;
	std::normal_distribution<double> hetDistro = std::normal_distribution<double>(mean,stddev);
	for(int i=0; i<this->nSystems; i++){
		double hetValue = -1.0;						// starting value is outside of acceptable range
		while(hetValue < 0.0 or hetValue < mean - 2.0*stddev or hetValue > mean + 2.0*stddev) hetValue = hetDistro(this->generator);		// loop until you find an acceptable value
		this->het[i] = hetValue;
	}
}



void Params::write_heterogeneity_values(){
	char buffer[300];
	std::ofstream dataout;
	sprintf(buffer,"u0_heterogeneity.bin");
	dataout.open(this->pthout+buffer,std::ios::binary);
	std::cout << "heterogeneity file: " << this->pthout+buffer << std::endl;
	if(!dataout){ std::cerr << red << "Error: Can not open output file " << this->pthout+buffer << reset << std::endl; exit(EXIT_FAILURE); }		// error checking
	for(int i=0; i<this->nSystems; i++) dataout.write((char*) &het[i], sizeof(double));
	dataout.close();
}




/// @brief constructor: set parameters
/// @param argc number of commandline arguments
/// @param argv values of commandline arguments
Params::Params(int argc, char *argv[]){
	
	// set up random numbers
	this->generator.seed(0);														// initialize generator
	this->normal_distro = std::normal_distribution<double>(0.0,1.0);				// random normal distribution with zero mean and unit standard deviation
	this->uniform_distro = std::uniform_real_distribution<double>(0.0,1.0);			// random uniform distribution in interval [0,1]
	
	// read parameterrs from commandline
	readCmdline(argc,argv);
	
	// parameters specified via commandline
	const double u0Dim = this->cmdParams[0];					// particle speed in m / s (20 micron/s)
	this->diffRot = this->cmdParams[1];							// diffusional rotation in rad^2 / s (see ...)
	this->diffS = this->cmdParams[2];							// noise intensity of run-tumble state (see ...)
	this->timeFactor = this->cmdParams[3];						// set times between jumps (match up)
	const double rDropletDim = this->cmdParams[4];				// droplet radius in m
	this->c = this->cmdParams[5];								// double well potential asymmetry factor
	this->swimFactor = this->cmdParams[6];						// swimming coefficient
	const double tfinal = this->cmdParams[7];					// total simulation time in s
	this->nSystems = this->cmdParams[8];						// number of systems with different random realizations
	this->noiseSeed = this->cmdParams[9];						// seed for deterministric noise generator
	this->hetQ = this->cmdParams[10];							// whether to use heterogeneity
	this->rateRT = this->cmdParams[11];							// R -> T rate in 1/s
	this->rateTR = this->cmdParams[12];							// T -> R rate in 1/s
	this->dragFudgeFactor = this->cmdParams[13];				// dragFudgeFactor
	
	// numerical parameters
	// OU
	// this->dt = 1.0e-2;											// = 1 s / 100 steps
	// this->nSaveStates = 1;										// how often to save state
	// RT
	this->dt = 1.0e-3;												// = 1 s / 1000 steps
	this->nSaveStates = 10;										// how often to save state
	// this->nSaveStates = 10;										// how often to save state
	// this->nSaveStates = 33;										// how often to save state -> output is generated at temporal resolution of 0.033s, as in experiments
	this->geometryChoice = model::confined_3d_torqueNoise_telegraph;	// good for RT on sphere								// 0 - sphere, 1 - 3d with projection, 2 - 3d + torqueNoise with projection, 3 - only bistable potential, 4 - 3d projection with telegraph
	// this->geometryChoice = model::planar_1d_ornstein_uhlenbeck;		// agrees with theory							// 0 - sphere, 1 - 3d with projection, 2 - 3d + torqueNoise with projection, 3 - only bistable potential, 4 - 3d projection with telegraph, 5 - 2d RT
	
	this->pthout = get_executable_path();								// path to executable = output directory
	std::cout << "pthout " << this->pthout << std::endl;
	this->nStates = tfinal/this->dt;
	this->nHistoryStates = this->nStates/this->nSaveStates + 1;			// account for initial condition via +1
    this->k = 2.0;     													// value of spring constant (Ornstein Uhlenbeck process)
	
	// physical parameters
	constexpr double gDim = 9.81;										// gravitational acceleration in m / s^2
	constexpr double rhoLCDim = 1.022e3;								// LC density in kg/m^3 (see paper p. 6)
	constexpr double rhoFCDim = 1.428e3;								// HFE 7200 density in kg/m^3 (see paper p. 6)
	constexpr double vR = 1.0;											// volume ratio of both phases: vR = vF / vH
	const double rInterfaceDim = 20.0*rDropletDim;						// inverse interface curvature in m
	constexpr double muDim = 8.9e-4;									// dynamic viscosity of water in Pa*s = kg / (m s)
	constexpr double rEcoliDim = 0.5e-6;								// short radius of E. coli bacteria in m
	constexpr double dragFactorT = 6.0*M_PI*muDim*rEcoliDim;			// translational Stokes drag factor of bacteria in kg/s; should be 8.38805e-9
	const double activeFrictionT = this->swimFactor * dragFactorT;		// active friction for swimming force
	// const double dragFactorR = 8.0*M_PI*muDim*pow(rDropletDim,3)/pow(rDropletDim+rEcoliDim,2);			// rotational Stokes drag factor of droplet in kg/s for force, see SI, not in torque, but linear force description
	const double dragFactorR = 8.0*M_PI*muDim*rDropletDim;				// rotational Stokes drag factor of droplet in kg/s for force, see SI, not in torque, but linear force description
	// const double dragFudgeFactor = 0.25;
	// const double dragFudgeFactor = 0.194545;	// average of ODE fit for speeds u = 4:2:20
	// const double dragFudgeFactor = 0.6;	// from literature
	// const double dragFudgeFactor = 0.5;
	// const double dragFudgeFactor = 1.0;
	const double dragFactor = dragFudgeFactor*(dragFactorT + dragFactorR);				// total drag due to bacteria + droplet
	
	// compute physical parameters
	const double mDim = calc_mass(rDropletDim,rhoLCDim,rhoFCDim,vR);							// droplet mass in kg
	const double dDim = calc_d(rDropletDim,rInterfaceDim,vR);									// distance d in m
	const double rcomDim = calc_rCOM(rDropletDim,rInterfaceDim,dDim,rhoLCDim,rhoFCDim,mDim);	// center of mass offset from geometrical center of sphere
	
	// std::cout << "dragFactorT: " << dragFactorT << " | 1.67761e-8 (rd=5e-6)" << std::endl;
	// std::cout << "dragFactorR: " << dragFactorR << " | ???" << std::endl;
	// std::cout << rDropletDim << " " << rInterfaceDim << " " << vR << std::endl;
	
	// rescaled parameters
	const double scale = rDropletDim;						// physical scaling	(for other scalings, r needs to be included explicitly in metric)
	this->rDroplet = rDropletDim / scale;					// rescaled droplet radius
	this->rForce = (rDropletDim+rEcoliDim) / scale;			// rescaled distance for droplet restoring force
	const double u0 = u0Dim / scale;						// rescaled particle speed
	const double g = gDim / scale;							// rescaled gravitational acceleration
	const double rcom = rcomDim / scale;					// rescaled center of mass offset
	
	// precompute parameters
	const double fudgeFactor = 1.0;
	this->swimFactor = activeFrictionT*u0/dragFactor;							// prefactor for swim force
	this->gravFactor = fudgeFactor*mDim*g*rcom/(this->rForce*dragFactor);		// prefactor for gravity force
	this->lcFactor = this->gravFactor;
	std::cout << this->swimFactor << " | " << this->gravFactor << std::endl;	// DEBUG
	
	// heterogeneity
	this->het.resize(this->nSystems);								// allocate heterogeneity array
	if(this->hetQ){													// fill heterogeneity factor array with heterogeneous values
		const double scalingFactor = activeFrictionT/(dragFactor*rDropletDim);	// swimFactor (=het) = u0 aFT / (drag rDropletDim)
		create_heterogeneity_values(scalingFactor);					// create heterogeneity values
		write_heterogeneity_values();								// write heterogeneity values to file
	}else{															// fill heterogeneity factor array with 1.0
		std::fill(het.begin(), het.end(), this->swimFactor);
	}
	
	// initialize number of variables
	if(int(this->geometryChoice) <= 7){
		const std::vector<int> nVariablesModels = {5,6,7,1,7,4,2,1};	// 
		this->nVariables = nVariablesModels[int(this->geometryChoice)];
	}else{
		printf("Error: Unknown geometry choice!\n"); 
		exit(EXIT_FAILURE);
	}
}

