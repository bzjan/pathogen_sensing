# pragma once

// functions for calculating droplet properties

double calc_mass(const double rDroplet, const double rhoF, const double rhoH, const double vR);
double calc_d(const double rd, const double ri, const double vr);
double calc_rCOM(const double rd, const double ri, const double d, const double rhoF, const double rhoH, const double m);

