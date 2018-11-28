/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef SPACE4ERRBODY_OPTIMIZATION_HPP
#define SPACE4ERRBODY_OPTIMIZATION_HPP

#include <utility>
#include <vector>
#include <limits>


#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/serialization.hpp"
#include "pagmo/problem.hpp"
#include "pagmo/types.hpp"

using namespace pagmo;
//using namespace pagmo::problem;

//! Define the problem PaGMO-style
struct Space4ErrBody {

    //! Empty constructor
    Space4ErrBody( ) { }

    //! Constructor.
    Space4ErrBody( const std::vector< std::vector< double > > &bounds,
                   const std::string &problem_name,
                   const std::string &vehicle_name,
                   const std::vector< std::string > &parameterList,
                   const std::vector< double > &parameterBounds,
                   const std::vector< double > &vehicleParameterValues,
                   const std::vector< std::string > &aeroCoeffFileList,
                   const std::vector< double > &simulation_settingsValues,
                   const std::vector< double > &initialConditionsValues,
                   const std::vector< double > &terminationConditionsValues,
                   const std::vector< double > &output_settingsValues,
                   const std::string &outputSubFolder );


    //! Calculate the fitness as a function of the parameter vector input_data
    std::vector< double > fitness( const std::vector< double >  &x ) const;

    //! Retrieve the allowable limits of the parameters to vary: pair
    //! containing minima and maxima of parameter values
    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    //! Retrieve the name of the problem
    std::string get_name( ) const;

    std::vector< double >::size_type get_nobj() const
    {
        return 9u;
    }
    std::vector< double >::size_type get_nec() const
    {
        return 0u;
    }

    //! Serialization function for Pagmo compatibility
    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

private:

    const std::vector< std::vector< double > > problemBounds_;
    const std::string problem_name_;
    const std::string vehicle_name_;
    const std::vector< std::string > parameterList_;
    const std::vector< double > parameterBounds_;
    const std::vector< double > vehicleParameterValues_;
    const std::vector< std::string > aeroCoeffFileList_;
    const std::vector< double > simulation_settingsValues_;
    const std::vector< double > initialConditionsValues_;
    const std::vector< double > terminationConditionsValues_;
    const std::vector< double >  output_settingsValues_;
    const std::string outputSubFolder_;

};


/*
namespace bislip {

struct parameters {

    Eigen::VectorXd xn_;
   std::string parameterName_;
    Eigen::VectorXd parameterValue_;

    //! Empty constructor
    parameters( ) {}

    //! Constructor
    parameters(
            Eigen::VectorXd xn,
            std::string const &parameterName,
            Eigen::VectorXd &parameterValue ):
        xn_(xn),
        parameterName_(parameterName),
        parameterValue_(parameterValue) {}

    bool operator < (const parameters &iData) const
    {
      return xn_ < iData.xn_;
    }



};

}*/

//PAGMO_REGISTER_PROBLEM(pagmo::Space4ErrBody)

#endif // SPACE4ERRBODY_OPTIMIZATION_HPP
