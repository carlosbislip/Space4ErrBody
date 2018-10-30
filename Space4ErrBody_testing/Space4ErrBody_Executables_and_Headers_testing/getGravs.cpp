#include "getStuff.h"

Eigen::VectorXd getGravs (
        const double mu,
        const double J2,
        const double J3,
        const double J4,
        const double R_E,
        const double r_norm,
        const double delta_rad ) {

    const double g_n_pre = -3 * ( mu / ( r_norm * r_norm ) ) * pow( R_E / r_norm , 2) * std::sin( delta_rad) * std::cos( delta_rad );
    const double g_n_sub_1 = J2;
    const double g_n_sub_2 = ( J3 / 2 ) * ( R_E / r_norm ) * ( 1 / std::sin( delta_rad) ) * ( 5 * pow (std::sin( delta_rad ) , 2 ) - 1 );
    const double g_n_sub_3 = ( 5 * J4 / 6 ) * pow( R_E / r_norm , 2) * ( 7 * pow (std::sin( delta_rad ) , 2 ) - 3 );
    const double g_n = g_n_pre * ( g_n_sub_1 + g_n_sub_2 + g_n_sub_3 );

    const double g_d_pre = ( mu / ( r_norm * r_norm ) );
    const double g_d_sub_1 = -( 3 * J2 / 2 ) * pow( R_E / r_norm , 2) * ( 3 * pow (std::sin( delta_rad ) , 2 ) - 1 );
    const double g_d_sub_2 = -2 * J3  * pow( R_E / r_norm , 3) * ( std::sin( delta_rad ) ) * ( 5 * pow (std::sin( delta_rad ) , 2 ) - 3 );
    const double g_d_sub_3 = -( 5 * J4 / 8 ) * pow( R_E / r_norm , 4 ) * ( 35 * pow (std::sin( delta_rad ) , 4 ) - 30 * pow (std::sin( delta_rad ) , 2 ) + 3 );
    const double g_d = g_d_pre * ( 1 + g_d_sub_1 + g_d_sub_2 + g_d_sub_3 );

    Eigen::VectorXd gravs ( 2 );
    gravs << g_n, g_d;

return gravs;

}
