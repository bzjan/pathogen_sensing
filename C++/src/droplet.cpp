// functions to compute droplet parameters

// C++ STL
#include <complex>								// std::complex<>
#include <vector>								// std::vector

// Eigen
#include <Eigen/Core>
#include <Eigen/Eigenvalues>					// Eigen::MatrixXd::eigenvalues()



/// @brief Calculate droplet mass m; [m] = kg. Source: Sara's thesis, Eq. (A.17) 
/// @param rDroplet droplet radius; [rDroplet] = m
/// @param rhoF density of fluorocarbon phase; [rhoF] = kg/m^3
/// @param rhoH density of hydrocarbon phase; [rhoH] = kg/m^3
/// @param vR volume ratio vR = vF/vH; [vR] = 1
/// @return droplet mass; [m] = kg
double calc_mass(const double rDroplet, const double rhoF, const double rhoH, const double vR){
	return 4.0/3.0*(rhoF + vR*rhoH)/(1.0+vR)*M_PI*pow(rDroplet,3);
}


/// @brief Finds roots of a poylnomial
/// @param coeffs Polynomial coefficients ordered from lowest order to highest, e.g.: coeffs = (a,b,c,d,e) with a x^0 + b x^1 + c x^2 + d x^3 + e x^4
/// @return vector containing roots of polynomial
std::vector<std::complex<double>> polynomial_root_finder(std::vector<double>& coeffs){
	
	const unsigned int deg = coeffs.size() - 1;
	std::vector<std::complex<double>> roots;
	Eigen::MatrixXd companion_mat = Eigen::MatrixXd::Zero(deg,deg);
	
	for(unsigned int n=0; n<deg; n++){
	for(unsigned int m=0; m<deg; m++){
		if(n == m + 1) companion_mat(n,m) = 1.0;
		if(m == deg-1) companion_mat(n,m) = -coeffs[n] / coeffs.back();
	}}
	
	Eigen::VectorXcd eigenValues = companion_mat.eigenvalues();
	for(unsigned int i=0; i<deg; i++) roots.push_back( eigenValues(i) );
	
	return roots;
}



/// @brief calculate center-center distance d
/// @param rd droplet radius
/// @param ri interface radius
/// @param vr volume ratio
/// @return distance d
double calc_d(const double rd, const double ri, const double vr){
	
	// get potential values for d (roots of 4th order polynomial)
	std::vector<double> coeffs = {-3.0*(pow(ri,4)+pow(rd,4)) + 6.0*pow(rd,2)*pow(ri,2), 8.0*(pow(ri,3)-pow(rd,3)) + 16.0*pow(rd,3)/(1.0+vr), -6.0*(pow(rd,2)+pow(ri,2)), 0.0, 1.0 };	// (a,b,c,d,e).(x^0,x^1,x^2,x^3,x^4)
	std::vector<std::complex<double>> roots = polynomial_root_finder(coeffs);
	
	// filter roots
	std::vector<double> filtered_roots;
	for(unsigned int i=0; i<roots.size(); i++){
		if(roots[i].imag() == 0.0 and -rd+ri < roots[i].real() and roots[i].real() < rd+ri ){
			filtered_roots.push_back( roots[i].real() );
		}
	}
	
	return filtered_roots[0];       // return first of found roots
}


/** Calculate center of mass distance to geometrical center of sphere. 
 * Works for any droplet geometry
 * @param rd droplet radius in m
 * @param ri interface radius in m
 * @param d center-center distance in m
 * @param rhoF density of fluorocarbon phase in kg/m^3
 * @param rhoH density of hydrocarbon phase in kg/m^3
 * @param m total mass of droplet in kg
 * @returns center of mass distance rCOM in m
 */
double calc_rCOM(const double rd, const double ri, const double d, const double rhoF, const double rhoH, const double m){
	const double ii = (pow(rd,2)-pow(ri,2)+pow(d,2))/(2.0*d);
	return (rhoH-rhoF)/(12.0*m)*M_PI*( pow(d-ri,2) * (d*d-3.0*ri*ri + 2.0*d*ri) + 3.0*pow(rd,4) - 4.0*d*pow(ii,3) );
}