// droplet function tests

// C++ STL
#include <cmath>                                // M_PI
#include <complex>								// std::complex<>
#include <vector>								// std::vector

// Eigen
#include <Eigen/Core>
#include <Eigen/Eigenvalues>					// Eigen::MatrixXd::eigenvalues()

// googletest
#include <gtest/gtest.h>                        // google test library
#include <gmock/gmock.h>                        // google mock library    

// custom 
#include "../src/droplet.h"                     // droplet functions


// use certain gmock functions
using ::testing::Pointwise, ::testing::FloatNear;


TEST(calc_mass, zeroResult){
    // test function output against known value
    constexpr double rDroplet = 0.0e-6;                     // droplet radius in m
    constexpr double rhoF = 1.01e3;									        // density in kg/m^3
    constexpr double rhoH = 1.44e3;									        // density in kg/m^3
    constexpr double vR = 1.0;											        // volume ratio of both phases: vR = vF / vH
    const double mCalc = calc_mass(rDroplet,rhoF,rhoH,vR);  // mass in kg
    EXPECT_NEAR(mCalc,0.0,1e-26);                           // compare values
}

TEST(calc_mass, correctSimpleResult) {
    // test function output against known value
    constexpr double rDroplet = 1.0;                     // droplet radius in m
    constexpr double rhoF = 1.0;									        // density in kg/m^3
    constexpr double rhoH = 1.0;									        // density in kg/m^3
    constexpr double vR = 1.0;											        // volume ratio of both phases: vR = vF / vH
    const double mCalc = calc_mass(rDroplet,rhoF,rhoH,vR);  // mass in kg
    const double mTrue = 4.0/3.0*M_PI;
    EXPECT_NEAR(mCalc,mTrue,1e-26);                           // compare values
}

TEST(calc_mass, correctPhysicalResult) {
    // test function output against known value
    constexpr double rDroplet = 3.5e-6;                     // droplet radius in m
    constexpr double rhoF = 1.01e3;									        // density in kg/m^3
    constexpr double rhoH = 1.44e3;									        // density in kg/m^3
    constexpr double vR = 1.0;											        // volume ratio of both phases: vR = vF / vH
    const double mCalc = calc_mass(rDroplet,rhoF,rhoH,vR);  // mass in kg
    constexpr double mMathematica = 2.200031155370152e-13;  // theoretically correct value
    EXPECT_NEAR(mCalc,mMathematica,1e-26);                  // compare values
}



// root finder

extern std::vector<std::complex<double>> polynomial_root_finder(std::vector<double>& coeffs);


TEST(polynomial_root_finder, correctValues) {
    
    // test true roots of poynomials
    std::vector<double> coeffs;
    std::vector<std::complex<double>> roots;
    std::vector<double> realRoots, trueRealRoots;
    auto const tolerance = 1e-5;
    
    // linear
	coeffs = { 0.0, 1.0 };	// f(x) = (x-0) = x
    trueRealRoots = {0.0};
	roots = polynomial_root_finder(coeffs);
    realRoots.resize(roots.size());
    for(unsigned int i=0; i<roots.size(); i++) realRoots[i] = roots[i].real();
    ASSERT_THAT(realRoots,Pointwise(FloatNear(tolerance), trueRealRoots));
    
    // quadratic
	coeffs = {-1.0, 0.0, 1.0 };	// f(x) = (x+1)(x-1) = x^2 - 1
    trueRealRoots = {-1.0,1.0};
	roots = polynomial_root_finder(coeffs);
    realRoots.resize(roots.size());
    for(unsigned int i=0; i<roots.size(); i++) realRoots[i] = roots[i].real();
    std::sort(realRoots.begin(), realRoots.end());
    ASSERT_THAT(realRoots,Pointwise(FloatNear(tolerance), trueRealRoots));
    
    // cubic
	coeffs = { 0.0, -1.0, 0.0, 1.0 };	// f(x) = (x+1)(x-0)(x-1) = x^3 - x
    trueRealRoots = {-1.0,0.0,1.0};
	roots = polynomial_root_finder(coeffs);
    realRoots.resize(roots.size());
    for(unsigned int i=0; i<roots.size(); i++) realRoots[i] = roots[i].real();
    std::sort(realRoots.begin(), realRoots.end());
    ASSERT_THAT(realRoots,Pointwise(FloatNear(tolerance), trueRealRoots));
    
    // quartic
	coeffs = { 0.0, -2.0, -1.0, 2.0, 1.0 };	// f(x) = (x-1)(x-0)(x+1)(x+2) = -2 x - x^2 + 2 x^3 + x^4
	trueRealRoots = {-2.0,-1.0,0.0,1.0};
    roots = polynomial_root_finder(coeffs);
    realRoots.resize(roots.size());
    for(unsigned int i=0; i<roots.size(); i++) realRoots[i] = roots[i].real();
    std::sort(realRoots.begin(), realRoots.end());
    ASSERT_THAT(realRoots,Pointwise(FloatNear(tolerance), trueRealRoots));
}




TEST(calc_d,correctPhysicalResult){
    // test function output against known value
    constexpr double rDroplet = 3.5e-6;                     // droplet radius in m
    constexpr double rInterface = 100.0e-6;					// density in kg/m^3
    constexpr double vR = 1.0;								// volume ratio of both phases: vR = vF / vH
    const double dCalc = calc_d(rDroplet,rInterface,vR);    // center-center distance
    constexpr double dMathematica = 0.0000999694;           // theoretically correct value
    EXPECT_NEAR(dCalc,dMathematica,1e-7);                   // compare values
}


TEST(calc_d,correctSimpleResult){
    // test function output against known value
    constexpr double rDroplet = 1.0;                        // droplet radius in m
    constexpr double rInterface = 100.0;					// density in kg/m^3
    constexpr double vR = 1.0;								// volume ratio of both phases: vR = vF / vH
    const double dCalc = calc_d(rDroplet,rInterface,vR);    // center-center distance
    constexpr double dMathematica = 99.9975;                // theoretically correct value
    EXPECT_NEAR(dCalc,dMathematica,1e-7);                   // compare values
}


TEST(calc_rCOM,correctPhysicalResult){
    // test function output against known value
    constexpr double rDroplet = 3.5e-6;                     // droplet radius in m
    constexpr double rInterface = 100.0e-6;					// inverse curvature of interface
    constexpr double vR = 1.0;								// volume ratio of both phases: vR = vF / vH
    constexpr double rhoF = 1.01e3;							// density in kg/m^3
    constexpr double rhoH = 1.44e3;							// density in kg/m^3
    const double m = calc_mass(rDroplet,rhoF,rhoH,vR);      // mass in kg
    const double d = calc_d(rDroplet,rInterface,vR);        // center-center distance
    const double rCOMCalc = calc_rCOM(rDroplet,rInterface,d,rhoF,rhoH,m);    // center-center distance
    constexpr double rCOMMathematica = 2.3034537919005271e-7;     // theoretically correct value
    ASSERT_LT(rCOMCalc,rDroplet) << "rCOM < rDroplet";
    ASSERT_NEAR(rCOMCalc,rCOMMathematica,1e-12);                   // compare values
    
    // test against vanishing curvature case (fluid phases are perfectly symmetric around center z=0 or theta=pi/2)
    const double rCOMTwoHalves = 3.0/8.0*rDroplet*(rhoH-rhoF)/(rhoH+rhoF);
    ASSERT_NEAR(rCOMCalc,rCOMTwoHalves,1e-7);                   // compare values
}