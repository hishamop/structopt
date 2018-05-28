#ifndef CMATERIAL_H
#define CMATERIAL_H
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"

namespace Material {
enum material_type{STEEL};

}

class CMaterial
{
public:
    CMaterial();
    ~CMaterial(){}
    virtual Eigen::Vector3d get_constitutive_function(Eigen::Vector3d,Eigen::Vector3d)=0;
    virtual const Eigen::Matrix3d& get_cmat()=0;
};

class STEEL:public CMaterial
{
public:
    STEEL();
    ~STEEL();
    Eigen::Vector3d get_constitutive_function(Eigen::Vector3d stress,Eigen::Vector3d strain){ return stress-m_stiffness_mat*strain;}
    const Eigen::Matrix3d& get_cmat(){return m_stiffness_mat;}

private:
    Eigen::Matrix3d m_stiffness_mat;
    double m_poisson;
    double m_young_modulus;
};



#endif // CMATERIAL_H
