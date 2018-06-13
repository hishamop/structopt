#include "cmaterial.h"

CMaterial::CMaterial()
{

}

STEEL::STEEL()
{
    m_poisson=0.3;
    m_young_modulus=2e5;
    double n=m_poisson;
    double E=m_young_modulus;
    Eigen::Matrix3d mat;
//    mat<<   1-n,n,0,
//            n,1-n,0,
//            0,0,1-2*n;
    mat<< 1/E,-n/E,0,
            -n/E,1/E,0,
            0,0,2*(1+n)/E;
    m_stiffness_mat= mat;
//    m_stiffness_mat= (m_young_modulus/((1+n)*(1-2*n)))*mat;
}
Eigen::Vector3d STEEL::get_constitutive_function(Eigen::Vector3d stress,Eigen::Vector3d strain)
{
//    return stress-m_stiffness_mat*strain;
    return strain- m_stiffness_mat*stress;
}

