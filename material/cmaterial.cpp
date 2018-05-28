#include "cmaterial.h"

CMaterial::CMaterial()
{

}

STEEL::STEEL()
{
    double n=m_poisson;
    Eigen::Matrix3d mat;
    mat<<   1-n,n,0,
            n,1-n,0,
            0,0,1-2*n;
    m_stiffness_mat= (m_young_modulus/((1+n)*(1-2*n)))*mat;
}
