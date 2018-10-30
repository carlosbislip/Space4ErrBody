#ifndef GET_STUFF_H
#define GET_STUFF_H

#include <iostream>
#include <ctime>
#include <sstream>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

std::vector< std::string > getFiles ( const std::string &filename );

std::vector< double > getData ( const std::string &filename );

std::string getCurrentDateTime ( const bool useLocalTime = false );

double getAngularDistance (
        const double &lat_i,
        const double &lon_i,
        const double &lat_f,
        const double &lon_f );

double getHeadingToTarget (
        const double &lat_i,
        const double &lon_i,
        const double &lat_f,
        const double &lon_f );


double getGamma ( const double A, const double B, const double phi );
double getEps_T ( const double &AoA, const Eigen::VectorXd &p, const Eigen::Vector6d &newCoefficients );


Eigen::VectorXd getGravs (
        const double mu,
        const double J2,
        const double J3,
        const double J4,
        const double R_E,
        const double r_norm,
        const double delta_rad );

double getRoot_NewtRaph (
        const double &x_old,
        const double &f,
        const double &fp);

double getRoot_Halley (
        const double &x_old,
        const double &f,
        const double &fp,
        const double &f2p);

double getRoot_Schroder (
        const double &x_old,
        const double &f,
        const double &fp,
        const double &f2p);

#endif // GET_STUFF_H
