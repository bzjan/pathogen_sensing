/*
 * Jan F. Totz (2021-2023)
 * program to simulate run-and-tumble bacteria attached to spherical droplets
 * 
 * runs on Windows & Ubuntu
 * 
 * compilation requires
 * * Eigen (>= 3.4) (https://eigen.tuxfamily.org/index.php)
 * 
 * cmake:
 * mkdir build
 * cd !$
 * cmake
 * 
 * 
 * Debug:
 * debug on linux: g++ -O0 -g -Wall -I /usr/include/eigen3/ run_tumble_on_sphere.cpp -o run_tumble_on_sphere.dbg.exe -fopenmp
 * valgrind --leak-check=full --track-origins=yes ./run_tumble_on_sphere.dbg.exe 2> valgrindoutput.txt
 */


// C++ STL
#include <vector>								// std::vector
#include <iostream>								// std::cout, std::cerr, std::endl
#include <fstream>								// std::ofstream
#include <random>								// default_random_engine, uniform_real_distribution
#include <algorithm>							// std::min



// custom
#include "./Params.h"							// Params
#include "./Timer.h"							// Timer
#include "./utility_functions.h" 				// red,blue,reset













// norm of 2-component vector on planar surface
double euclidean_norm(const double v1, const double v2){
	return sqrt(v1*v1 + v2*v2);
}


// norm of 3-component vector on planar surface
double euclidean_norm(const double v1, const double v2, const double v3){
	return sqrt(v1*v1 + v2*v2 + v3*v3);
}


// norm on tangent space of unit sphere
// metric g = {{1,0},{0,sin^2 theta}}
double spherical_norm(const double v1, const double v2, const double theta){
	return sqrt(v1*v1 + sin(theta)*sin(theta) * v2*v2);
}


void apply_euclidean_norm_constraint(double& pTheta, double& pPhi){
	const double norm = euclidean_norm(pTheta,pPhi);
	pTheta /= norm;
	pPhi /= norm;
}


void apply_euclidean_norm_constraint(double& x, double& y, double& z){
	const double norm = euclidean_norm(x,y,z);
	x /= norm;
	y /= norm;
	z /= norm;
}

void apply_euclidean_2d_norm_constraint(double& x, double& y){
	const double norm = euclidean_norm(x,y);
	x /= norm;
	y /= norm;
}

// TODO: improve! right now incorrect! -> weight components differently!
void apply_spherical_norm_constraint(double& pTheta, double& pPhi, const double theta){
	const double norm = spherical_norm(pTheta,pPhi,theta);
	pTheta /= norm;
	pPhi /= norm;
}


void save_state(const std::vector<double>& vars, std::vector<double>& histories, int& saveCounter, const Params& p){
	
	const int tOffset = saveCounter*p.nVariables*p.nSystems;
	for(int i=0; i<p.nSystems; i++){
	for(int j=0; j<p.nVariables; j++){
		histories[j + i*p.nVariables + tOffset] = vars[i+j*p.nSystems];
	}}
	saveCounter++;
}

double pow3(const double x){
	return x*x*x;
}

// derivative of bistable potential
double duds(const double x, const double c){
	// double potential parameters
	const double a = 32.0;
	const double b = 1.0;
	return 11.3*c - 11.3*a*(-6.0 + 11.3*x) + 11.3*b*pow3(-6.0 + 11.3*x);
}


// TODO: buggy; redo! Castro-Villarreal 2018
void dynamics_sphere(double& theta, double& phi, double& pTheta, double& pPhi, double& s, Params& p){
	
	const double noiseTheta = p.normal_distro(p.generator);
	const double noisePhi = p.normal_distro(p.generator);
	const double noiseS = p.normal_distro(p.generator);
	
	// dynamic evolution
	const double phiDot = (1.0-s)*p.swimFactor*pPhi/(p.rForce*sin(theta));								// phi_t = phi + dt * phiDot, with phiDot = f_phi(theta,phi)
	const double thetaNew = theta + p.dt*( (1.0-s)*p.swimFactor*pTheta/p.rForce - p.gravFactor*sin(theta) );
	const double phiNew = phi + p.dt*phiDot;
	const double pThetaNew = pTheta + sqrt(p.dt)*( sqrt(2.0*s*p.diffRot)*pPhi*( pPhi*noiseTheta - pTheta*noisePhi ) + cos(theta)*pPhi*phiDot );		// theta component of polarization vector p = pTheta*eTheta + pPhi*ePhi
	const double pPhiNew = pPhi + sqrt(p.dt)*( sqrt(2.0*s*p.diffRot)*pTheta*( pTheta*noisePhi - pPhi*noiseTheta ) - cos(theta)*pTheta*phiDot );		// phi component of polarization vector p = pTheta*eTheta + pPhi*ePhi
	const double sNew = s + p.dt*( -p.timeFactor*duds(s,p.c) ) + sqrt(p.dt*p.timeFactor*p.diffS)*noiseS;																		// // run/tumble state variable
	
	// update variables
	theta = thetaNew;
	phi = phiNew;
	pTheta = pThetaNew;
	pPhi = pPhiNew;
	s = std::min(std::max(sNew,0.0),1.0);			// super important! otherwise negative diffusion -> instability

	// project back on sphere to avoid numerical issues
	apply_euclidean_norm_constraint(pTheta,pPhi);		// reset unit polarization vector length: ||p|| = 1
}


// continuous-time two state random process
// source: Lindner chapter in Laing, Lord 2009
void telegraph_process(double& s, Params& p){ 
	const double probability = p.uniform_distro(p.generator);
	if(s==0.0){		// state s=0 (running)
		if(probability < p.rateRT*p.dt) s = 1.0;
	}else{			// state s=1 (tumbling)
		if(probability < p.rateTR*p.dt) s = 0.0;
	}
}


// cross product: c = r x p
// output overwrites px,py,pz
void crossProduct(
const double x, const double y, const double z, 
const double px, const double py, const double pz, 
double& cx, double& cy, double& cz){
	cx = y*pz - z*py;
	cy = z*px - x*pz;
	cz = x*py - z*px;
}


/// @brief project force vector to local tangential plane on sphere
/// @param x position vector, x-component
/// @param y position vector, y-component
/// @param z position vector, z-component
/// @param vx orientation/director vector, x-component
/// @param vy orientation/director vector, y-component
/// @param vz orientation/director vector, z-component
void project_force_to_sphere(
	const double x, const double y, const double z, 
	double& vx, double& vy, double& vz
){
	// assume position is on sphere
	
	// 1) surface projection operator: P_S p = (1 - n \otimes n) p = p - (n.p) n with n = {x,y,z} unit surface normal vector (projected on sphere above)
	const double nDotp = x*vx + y*vy + z*vz;			// dot product
	vx -= nDotp*x;
	vy -= nDotp*y;
	vz -= nDotp*z;
}


/// @brief project gravity force vector to local tangential plane on sphere with delayed onset
/// @param theta position vector, theta coordinate
/// @param phi position vector, phi coordinate
/// @param gravFactor gravitational force prefactor
/// @param fx force vector, x-component
/// @param fy force vector, y-component
/// @param fz force vector, z-component
void project_force_to_sphere_grav(
	const double theta, const double phi,
	const double gravFactor,
	double& fx, double& fy, double& fz
){	
	// surface projection operator for force: F = P_S e_z with P_S = 1 - e_r \otimes e_r
	fx += -gravFactor * (cos(theta)*cos(phi)*sin(theta));
	fy += -gravFactor * (cos(theta)*sin(theta)*sin(phi));
	fz +=  gravFactor * sin(theta)*sin(theta);
}


/// @brief project free trajectory back to position sphere and orientation to orientation sphere
/// @param x position vector, x-component
/// @param y position vector, y-component
/// @param z position vector, z-component
/// @param px orientation/director vector, x-component
/// @param py orientation/director vector, y-component
/// @param pz orientation/director vector, z-component
void project_to_sphere(double& x, double& y, double& z, double& px, double& py, double& pz){
	
	// project position vector to sphere: r = r / |r|
	double norm = euclidean_norm(x,y,z);
	x /= norm;
	y /= norm;
	z /= norm;
	
	// 1) surface projection operator: P_S p = (1 - n \otimes n) p = p - (n.p) n with n = {x,y,z} unit surface normal vector (projected on sphere above)
	const double nDotp = x*px + y*py + z*pz;			// dot product
	px -= nDotp*x;
	py -= nDotp*y;
	pz -= nDotp*z;
	
	// normalize polarization vector: p = p / |p|
	norm = euclidean_norm(px,py,pz);
	px /= norm;
	py /= norm;
	pz /= norm;
}


void force_liquid_crystal(
	const double theta, const double phi,
	const double kLC,
	double& fx, double& fy, double& fz
){
	// F_LC = -kLC * (thetaLC-0) e_theta
	fx += -kLC * theta * (cos(theta)*cos(phi));
	fy += -kLC * theta * (cos(theta)*sin(phi));
	fz += -kLC * theta * (-sin(theta));
}



// calculations are in dimensionless coordinates
void dynamics_confined_3d_torqueNoise_telegraph(
	double& x, double& y, double& z, 
	double& px, double& py, double& pz, 
	double& s, 
	const double het, Params& p
){
	// dynamic evolution: dot{r} = P(F_swim + F_tau) = P(u0 p(pTheta,pPhi) + mDroplet*g*rCOM/rForce / dragFactor e_z)
	
	const double phi = atan2(y,x);
	const double theta = acos(z);						// z in dimensionless coordinates
	// constexpr double thetaThresh = 10.0*M_PI/180.0;		// from degrees to radians
	// constexpr double thetaThresh = 0.0*M_PI/180.0;		// from degrees to radians
	// const double thetaEffective = theta - thetaThresh;
	
	// right hand side of xdot = F (projected force components)
	double fx = 0.0;
	double fy = 0.0;
	double fz = 0.0;
	
	if(s==0.0){		// running s = 0
	
		// calculate projected swim force P(F_swim)
		fx = het*px;								// F_swim_x
		fy = het*py;								// F_swim_y
		fz = het*pz;								// F_swim_z
		project_force_to_sphere(x,y,z,fx,fy,fz);	// P(F_swim)
		
	}else{			// tumbling s = 1
		
		// dot{p} = T = noise*(r x p)
		// equivalent to rotating director around surface normal
		double cx,cy,cz;								// cross product components (torque T)
		crossProduct(x,y,z,px,py,pz,cx,cy,cz);
		const double noiseTorque = sqrt(p.dt*2*p.diffRot)*p.normal_distro(p.generator);		// diffusionCoeff*torque = rotationalDiffusion
		
		// update director components
		px += noiseTorque*cx;
		py += noiseTorque*cy;
		pz += noiseTorque*cz;
	}
	
	// calculate projected gravity force P(F_g)
	// if(thetaEffective < 0){
	// 	project_force_to_sphere_grav(theta,phi,p.gravFactor,fx,fy,fz);
	// }else{
	// 	project_force_to_sphere_grav(thetaEffective,phi,0.5*p.gravFactor,fx,fy,fz);
	// }
		
	project_force_to_sphere_grav(theta,phi,p.gravFactor,fx,fy,fz);
	// project_force_to_sphere_grav(0.5*theta,phi,p.gravFactor,fx,fy,fz);
	
	// calculate LC force "P(F_LC)"
	// force_liquid_crystal(0.5*theta,phi,p.lcFactor,fx,fy,fz);
	// if(thetaEffective > 0) force_liquid_crystal(0.5*thetaEffective,phi,p.gravFactor,fx,fy,fz);
	
	// update position variables: xnew = x + dt*F (overdamped motion)
	x += p.dt*fx;
	y += p.dt*fy;
	z += p.dt*fz;
	project_to_sphere(x,y,z,px,py,pz);		// project {position, direction} state onto (sphere x sphere)
	
	// update stochastic state via telegraph process for switching between run and tumble states
	telegraph_process(s,p);		// update s in place to avoid memory usage (no other variables couple into s)
}



// calculations are in dimensionless coordinates
void dynamics_spherical_gravity_liquid_crystal(
	double& x, double& y, double& z, 
	double& px, double& py, double& pz, 
	double& s, 
	double& thetaG, double& thetaLC,
	const double het, Params& p
){
	// dynamic evolution: dot{r} = P(F_swim + F_tau) = P(u0 p(pTheta,pPhi) + mDroplet*g*rCOM/rForce / dragFactor e_z)
	
	const double phi = atan2(y,x);
	// constexpr double thetaThresh = 10.0*M_PI/180.0;		// from degrees to radians
	constexpr double thetaGThresh = 10.0*M_PI/180.0;		// from degrees to radians
	const double thetaGEffective = thetaG - thetaGThresh;
	
	// right hand side of xdot = F (projected force components)
	double fx = 0.0;
	double fy = 0.0;
	double fz = 0.0;
	
	if(s==0.0){		// running s = 0
	
		// calculate projected swim force P(F_swim)
		fx = het*px;								// F_swim_x
		fy = het*py;								// F_swim_y
		fz = het*pz;								// F_swim_z
		project_force_to_sphere(x,y,z,fx,fy,fz);	// P(F_swim)
		
	}else{			// tumbling s = 1
		
		// dot{p} = T = noise*(r x p)
		// equivalent to rotating director around surface normal
		double cx,cy,cz;								// cross product components (torque T)
		crossProduct(x,y,z,px,py,pz,cx,cy,cz);
		const double noiseTorque = sqrt(p.dt*2*p.diffRot)*p.normal_distro(p.generator);		// diffusionCoeff*torque = rotationalDiffusion
		
		// update director components
		px += noiseTorque*cx;
		py += noiseTorque*cy;
		pz += noiseTorque*cz;
	}
	
	// calculate projected gravity force P(F_g)
	if(thetaGEffective > 0) project_force_to_sphere_grav(thetaGEffective,phi,p.gravFactor,fx,fy,fz);
	
	// calculate LC force "P(F_LC)"
	force_liquid_crystal(thetaLC,phi,p.lcFactor,fx,fy,fz);
	
	// update position variables: xnew = x + dt*F (overdamped motion)
	x += p.dt*fx;
	y += p.dt*fy;
	z += p.dt*fz;
	
	
	// update stochastic state via telegraph process for switching between run and tumble states
	telegraph_process(s,p);		// update s in place to avoid memory usage (no other variables couple into s)
	
	// project state onto sphere x sphere
	project_to_sphere(x,y,z,px,py,pz);
}





void dynamics_planar_2d_torqueNoise_telegraph(
	double& x, double& y,  
	double& omega, 
	double& s, 
	const double het, Params& p
){
	const double noiseOmega = p.normal_distro(p.generator);
	
	// dynamic evolution: dot{r} = F_swim + F_tau = u0 p(pTheta,pPhi) + mDroplet*g*rCOM/rForce / dragFactor e_z
	const double xNew = x + p.dt*( (1.0-s)*het*cos(omega) );
	const double yNew = y + p.dt*( (1.0-s)*het*sin(omega) );
	const double omegaNew = omega + sqrt(p.dt*2*p.diffRot*s)*noiseOmega;
	
	// telegraph process for switching between run and tumble states
	telegraph_process(s,p);		// update s in place to avoid memory usage (no other variables couple into s)
	
	// update variables
	x = xNew;
	y = yNew;
	omega = omegaNew;
}



/// @brief 1d Ornstein Uhlenbeck stochastic dynamics
/// @param x position
/// @param p parameters
void dynamics_1d_ornstein_uhlenbeck(
	double& x, 
	Params& p
){
	const double noiseX = p.normal_distro(p.generator);
	
	// dynamic evolution: dot{r} = -k*r + D xi + update
	x += p.dt*( -p.k*x ) + sqrt(p.dt)*noiseX;
}


/// @brief 1d Wiener process stochastic dynamics
/// @param x position
/// @param p parameters
void dynamics_1d_wiener(
	double& x, 
	Params& p
){
	// dynamic evolution: dot{r} = D xi
	x += sqrt(p.dt)*p.normal_distro(p.generator);
}


/// @brief 2d Ornstein Uhlenbeck stochastic dynamics
/// @param x position
/// @param y position
/// @param p parameters
void dynamics_2d_ornstein_uhlenbeck(
	double& x, double& y,  
	Params& p
){
	const double noiseX = p.normal_distro(p.generator);
	const double noiseY = p.normal_distro(p.generator);
	
	// dynamic evolution: dot{r} = -k*r + D xi + update
	x += p.dt*( -p.k*x ) + sqrt(p.dt)*noiseX;
	y += p.dt*( -p.k*y ) + sqrt(p.dt)*noiseY;
}


// projected 3d dynamics with gravity
void dynamics_confined_3d(double& x, double& y, double& z, double& pTheta, double& pPhi, double& s, Params& p){
	
	const double noise1 = p.normal_distro(p.generator);
	const double noise2 = p.normal_distro(p.generator);
	const double noise3 = p.normal_distro(p.generator);
	const double noiseS = p.normal_distro(p.generator);
	
	// dynamic evolution: dot{r} = F_swim + F_tau = u0 e_r(pTheta,pPhi) + mDroplet*g*rCOM/rForce / dragFactor e_z
	// Note: projection to sphere happens in a later step 
	const double xNew = x + p.dt*( (1.0-s)*p.swimFactor*cos(pPhi)*sin(pTheta) );
	const double yNew = y + p.dt*( (1.0-s)*p.swimFactor*sin(pPhi)*sin(pTheta) );
	const double zNew = z + p.dt*( (1.0-s)*p.swimFactor*cos(pTheta) + p.gravFactor );
	// brownian motion on sphere for pTheta and pPhi (polarization unit vector p = e_r(pTheta,pPhi))
	// pThetaDot = 0.5*Cot[pTheta]*diffR + sqrt(diffR)*(Sin[pPhi]*w1 - Cos[pPhi]*w2)
	// pPhiDot = sqrt(diffR)*(cot[pTheta]*(cos[pPhi]*w1 + sin[pPhi]*w2) - w3)
	const double pThetaNew = pTheta + p.dt*( s*p.diffRot*0.5*cos(pTheta)/sin(pTheta) ) + sqrt(p.dt)*( sqrt(s*p.diffRot)*(sin(pPhi)*noise1 - cos(pPhi)*noise2) );
	const double pPhiNew = pPhi + sqrt(p.dt)*( sqrt(s*p.diffRot)*( cos(pTheta)/sin(pTheta)*(cos(pPhi)*noise1 + sin(pPhi)*noise2) - noise3 ) );
	// double potential evolution for switching between run and tumble states
	const double sNew = s + p.dt*( -p.timeFactor*duds(s,p.c) ) + sqrt(p.dt*p.timeFactor*p.diffS)*noiseS;
	
	// update variables
	x = xNew;
	y = yNew;
	z = zNew;
	pTheta = pThetaNew;
	pPhi = pPhiNew;
	s = std::min(std::max(sNew,0.0),1.0);			// super important! otherwise negative diffusion -> instability
}




// confined 3d dynamics with gravity
// as in Sknepnek, Henkes 2015
void dynamics_confined_3d_torqueNoise(
	double& x, double& y, double& z, 
	double& px, double& py, double& pz, double& s,
	const double het, Params& p
){
	const double noiseTorque = p.normal_distro(p.generator);
	const double noiseS = p.normal_distro(p.generator);
	
	// dynamic evolution: dot{r} = F_swim + F_tau = u0 p(pTheta,pPhi) + mDroplet*g*rCOM/rForce / dragFactor e_z
	// Note: projection to sphere happens in a later step 
	// het = p.swimFactor
	const double xNew = x + p.dt*( (1.0-s)*het*px );
	const double yNew = y + p.dt*( (1.0-s)*het*py );
	const double zNew = z + p.dt*( (1.0-s)*het*pz + p.gravFactor );
	
	// dot{p} = noise*(r x p)
	// equivalent to rotating director around surface normal
	double cx,cy,cz;								// cross product components
	crossProduct(x,y,z,px,py,pz,cx,cy,cz);
	const double pxNew = px + sqrt(p.dt*2*p.diffRot*s)*noiseTorque*cx;
	const double pyNew = py + sqrt(p.dt*2*p.diffRot*s)*noiseTorque*cy;
	const double pzNew = pz + sqrt(p.dt*2*p.diffRot*s)*noiseTorque*cz;
	
	// double potential evolution for switching between run and tumble states
	const double sNew = s + p.dt*( -p.timeFactor*duds(s,p.c) ) + sqrt(p.dt*p.timeFactor*p.diffS)*noiseS;
	//~ const double sNew = s;
	
	// update variables
	x = xNew;
	y = yNew;
	z = zNew;
	px = pxNew;
	py = pyNew;
	pz = pzNew;
	s = std::min(std::max(sNew,0.0),1.0);			// super important! otherwise negative diffusion -> instability
}



/// @brief stochastic dynamics in bistable potential
/// @param s stochastic state
/// @param p parameter vector
void dynamics_bistable_potential(double& s, Params& p){
	
	const double noiseS = p.normal_distro(p.generator);
		
	// double potential evolution for switching between run and tumble states
	const double sNew = s + p.dt*( -p.timeFactor*duds(s,p.c) ) + sqrt(p.dt*p.timeFactor*p.diffS)*noiseS;
	
	// update variables
	s = std::min(std::max(sNew,0.0),1.0);			// super important! otherwise negative diffusion -> instability
}




/// @brief write output to file
/// @param histories 
/// @param p 
void write_output_to_file(const std::vector<double>& histories, const Params& p){
	
	std::ofstream dataout;
	const std::string pthfn = p.pthout+"runTumbleOutput.bin";	// pth / filename
	dataout.open(pthfn,std::ios::binary);
	std::cout << pthfn << std::endl;
	if(!dataout){ std::cerr << red << "Error: Can not open output file " << pthfn << reset << std::endl; exit(EXIT_FAILURE); }		// error checking
	
	for(int j=0; j<p.nSystems; j++){
	for(int t=0; t<p.nHistoryStates; t++){
	for(int i=0; i<p.nVariables; i++){
		dataout.write((char*) &histories[i + j*p.nVariables + t*p.nVariables*p.nSystems], sizeof(double));
	}}}
	
	dataout.close();
}











/// @brief returns modulo with positive result
/// @param x value
/// @param y divisor
/// @return x % y in R+
inline double posi_fmod(const double x, const double y){ return fmod((fmod(x,y)+y),y); }

// project free trajectory back to sphere
void project_to_sphere(double& x, double& y, double& z, double& pTheta, double& pPhi){
	
	// project position vector to sphere: r = r / |r|
	const double norm = euclidean_norm(x,y,z);
	x /= norm;
	y /= norm;
	z /= norm;
	
	// project polarization vector p = e_r(pTheta,pPhi) on tangent plane of current position on unit sphere
	double px = cos(pPhi)*sin(pTheta);
	double py = sin(pPhi)*sin(pTheta);
	double pz = cos(pTheta);
	
	// 1) surface projection operator: P_S p = (1 - n \otimes n) p = p - (n.p) n with n = {x,y,z} unit surface normal vector (projected on sphere above)
	const double nDotp = x*px + y*py + z*pz;			// dot product
	px -= nDotp*x;
	py -= nDotp*y;
	pz -= nDotp*z;
	
	// 2) exctract angles from projected vector (convert from Cartesian to spherical coordinates)
	pTheta = -atan2(pz,euclidean_norm(px,py)) + 0.5*M_PI;			// returns [0,pi] from north to south
	pPhi = posi_fmod(atan2(py,px),2.0*M_PI);						// returns [0,2pi] ccw from x to y axis
}



void geometry_confined_3d(std::vector<double>& vars, Params& p){
	
	// multiple random realizations
	// #pragma omp parallel for
	for(int i=0; i<p.nSystems; i++){
		dynamics_confined_3d(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems],p);	// evolve in time (overdamped particle)
		project_to_sphere(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems]);					// project state onto sphere
	}
}

void geometry_confined_3d_torqueNoise(std::vector<double>& vars, Params& p){
	
	// multiple random realizations
	// #pragma omp parallel for
	for(int i=0; i<p.nSystems; i++){
		dynamics_confined_3d_torqueNoise(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems],vars[i+6*p.nSystems],p.het[i],p);		// evolve in time (overdamped particle)
		project_to_sphere(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems]);									// project state onto sphere
	}
}

// double calc_thetaNoise(Params& p){

// 	return thetaDot + mg rRatio sin(theta);
// }



/// @brief set initial conditions
/// @param vars 
/// @param p 
void initial_conditions(std::vector<double>& vars, Params& p){
	
	std::vector<double> ic;
	std::default_random_engine g(p.noiseSeed);						// seed random number generator
	std::uniform_real_distribution<double> uniform_p(-1.0,1.0);		// uniform random numbers between -1 and 1 for polarization vector components
	std::uniform_real_distribution<double> uniform_s(0.0,1.0);		// uniform random numbers between 0 and 1 for random state
	
	// start all realizations from the same initial condition and then let them evolve differently
	switch(p.geometryChoice){
		case model::spherical:																// sphere
			for(int i=0; i<p.nSystems; i++){
				ic = {0.5*M_PI,0.1, uniform_s(g)*M_PI,uniform_s(g)*2.0*M_PI, uniform_s(g)};						// (spherical position:) theta,phi, (angular direction of polarization vector:) pTheta,pPhi, s
				if(ic[0]==0.0){ printf("Warning: Polar angle theta must NOT be zero! Correcting now to 0.5 pi\n"); ic[0]=0.5*M_PI; }
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
				apply_euclidean_norm_constraint(vars[i+2*p.nSystems],vars[i+3*p.nSystems]);
			}
			break;

		case model::confined_3d:																// confined 3d space
			for(int i=0; i<p.nSystems; i++){
				ic = { 0.0,1.0,1.0, uniform_s(g)*M_PI, uniform_s(g)*2.0*M_PI, uniform_s(g) };					// x,y,z, (angular direction of polarization vector:) pTheta,pPhi, s
				if(ic[3]==0.0){ printf("Warning: Polar angle pTheta must NOT be zero! Correcting now to 0.5 pi\n"); ic[3]=0.5*M_PI; }
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
				project_to_sphere(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems]);		// project state onto sphere
			}
			break;
			
		case model::confined_3d_torqueNoise:																// confined 3d space with torque
			for(int i=0; i<p.nSystems; i++){
				ic = {0.0,1.0,1.0, uniform_p(g),uniform_p(g),uniform_p(g), uniform_s(g)};				// x,y,z, (cartesian components of polarization vector:) px,py,pz, s
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
				project_to_sphere(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems]);		// project state onto sphere
			}
			break;
			
		case model::bistable_potential:																// only stochastic variable dynamics
			for(int i=0; i<p.nSystems; i++){
				vars[i] = uniform_s(g);										//  s
			}
			break;

		case model::confined_3d_torqueNoise_telegraph:																// confined 3d space with telegraph noise for RT
			for(int i=0; i<p.nSystems; i++){
				ic = {0.0,1.0,1.0, uniform_p(g),uniform_p(g),uniform_p(g), 0.0};				// x,y,z, (cartesian components of polarization vector:) px,py,pz, s
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
				project_to_sphere(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems]);		// project state onto sphere
			}
			break;
			
		case model::planar_2d_torqueNoise_telegraph:																// 2d RT motion
			for(int i=0; i<p.nSystems; i++){
				ic = {0.0,1.0, uniform_p(g), 0.0};				// x,y, (cartesian components of polarization vector:) px,py, s
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
				apply_euclidean_2d_norm_constraint(vars[i+2*p.nSystems],vars[i+3*p.nSystems]);		// normalize polarization vector
			}
			break;
			
		case model::planar_2d_ornstein_uhlenbeck:										// 2d OU
			for(int i=0; i<p.nSystems; i++){
				ic = {0.0,1.0};						// x,y
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
			}
			break;
			
		case model::planar_1d_ornstein_uhlenbeck:										// 1d OU
			for(int i=0; i<p.nSystems; i++){
				ic = {1.0};							// x
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
			}
			break;
		
		case model::spherical_gravity_liquid_crystal:										// 1d OU
			for(int i=0; i<p.nSystems; i++){
				ic = {0.0,0.0,1.0, uniform_p(g),uniform_p(g),uniform_p(g), uniform_s(g),0.0,0.0};							// x,y,z,px,py,pz,s,thetaG,thetaLC
				for(int j=0; j<p.nVariables; j++) vars[i+j*p.nSystems] = ic[j];
			}
			break;
					
		default:
			printf("Error: Unknown geometry choice!\n"); exit(EXIT_FAILURE);
			break;
	}
}

// time evolution
void timeEvolution(std::vector<double>& vars, std::vector<double>& histories, int& saveCounter, Params& p){
	
	for(int t=0; t<p.nStates; t++){
		switch(p.geometryChoice){
			case model::spherical:
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_sphere(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],p); 
				break;
			case model::confined_3d: geometry_confined_3d(vars,p); break;
			case model::confined_3d_torqueNoise: geometry_confined_3d_torqueNoise(vars,p); break;
			case model::bistable_potential:
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_bistable_potential(vars[i],p);
				break;
			case model::confined_3d_torqueNoise_telegraph: 	// good version!
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_confined_3d_torqueNoise_telegraph(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems],vars[i+6*p.nSystems],p.het[i],p);
				break;
			case model::planar_2d_torqueNoise_telegraph: 
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_planar_2d_torqueNoise_telegraph(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],p.het[i],p);
				break;
			case model::planar_2d_ornstein_uhlenbeck:
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_2d_ornstein_uhlenbeck(vars[i],vars[i+p.nSystems],p);
				break;
			case model::planar_1d_ornstein_uhlenbeck: 
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_1d_ornstein_uhlenbeck(vars[i],p);
				break;
			case model::spherical_gravity_liquid_crystal:
				// #pragma omp parallel for
				for(int i=0; i<p.nSystems; i++) dynamics_spherical_gravity_liquid_crystal(vars[i],vars[i+p.nSystems],vars[i+2*p.nSystems],vars[i+3*p.nSystems],vars[i+4*p.nSystems],vars[i+5*p.nSystems],vars[i+7*p.nSystems],vars[i+8*p.nSystems],vars[i+9*p.nSystems],vars[i+10*p.nSystems],p);
				break;
			default: printf("Error: Unknown dynamics choice!\n"); exit(EXIT_FAILURE); break;
		}
		// save state
		if( !(t%p.nSaveStates) ) save_state(vars,histories,saveCounter,p);
	}
}





/// @brief main function
/// @param argc number of input arguments
/// @param argv values of input arguments
/// @return program status
int main(int argc, char *argv[]){
	
	/// 1) set parameters
	Params p(argc,argv);														// parameter container
	
	/// 2) initialization
	
	// state vectors
	std::vector<double> vars(p.nSystems*p.nVariables);							// vars = {theta[All], phi[All], ...}
	std::vector<double> histories(p.nVariables*p.nSystems*p.nHistoryStates);	// history vector
	
	// initial condition
	initial_conditions(vars,p);
	
	// save initial state
	int saveCounter = 0;
	save_state(vars,histories,saveCounter,p);
	
	// timer
	Timer timer;
	timer.start();										// timer start
	
	/// 3) simulation
	timeEvolution(vars,histories,saveCounter,p);		// evolve temporal dynamics
	
	timer.end();										// timer end
	timer.print_runtime("Computation finished");		// terminal output
	
	
	/// 4) save data to file
	
	// write simulation
	write_output_to_file(histories,p);
	
	return EXIT_SUCCESS;
}

