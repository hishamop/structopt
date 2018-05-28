#include "shapefn.h"

void Quad4::update_shapefn( double sval, double tval)
{
  //assert(m_xcoord.size()==4 and m_ycoord.size()==4);
  m_shapefn_val.resize(4);
  m_Nis.resize(4);
  m_Nit.resize(4);

  m_shapefn_val(0)=0.25*(1-sval)*(1-tval);
  m_shapefn_val(1)=0.25*(1+sval)*(1-tval);
  m_shapefn_val(2)=0.25*(1+sval)*(1+tval);
  m_shapefn_val(3)=0.25*(1-sval)*(1+tval);

  m_Nis[0]=(tval-1)/4;
  m_Nis[1]=-(tval-1)/4;
  m_Nis[2]=(tval+1)/4;
  m_Nis[3]=-(tval+1)/4;

  m_Nit[0]= (sval-1)/4;
  m_Nit[1]=-(sval+1)/4;
  m_Nit[2]=(sval+1)/4;
  m_Nit[3]=-(sval-1)/4;

  Eigen::Map<Eigen::VectorXd>m_x(m_xcoord.data(),m_xcoord.size());
  Eigen::Map<Eigen::VectorXd>m_y(m_ycoord.data(),m_ycoord.size());
  m_jacobian.resize(2,2);
  m_jacobian(0,0) = m_Nis*m_x;
  m_jacobian(0,1) = m_Nit*m_x;
  m_jacobian(1,0) = m_Nis*m_y;
  m_jacobian(1,1) = m_Nit*m_y;

  m_inverse_jacobian.resize(2,2);
  m_detJ= m_jacobian(0,0)*m_jacobian(1,1)-m_jacobian(0,1)*m_jacobian(1,0);
  m_inverse_jacobian(0,0)=m_detJ*m_jacobian(1,1);
  m_inverse_jacobian(0,1)=-1*m_detJ*m_jacobian(0,1);
  m_inverse_jacobian(1,0)=-1*m_detJ*m_jacobian(1,0);
  m_inverse_jacobian(1,1)=m_detJ*m_jacobian(0,0);
}







