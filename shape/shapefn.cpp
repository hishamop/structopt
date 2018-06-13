#include "shapefn.h"

void Quad4::update_shapefn( double sval, double tval)
{
  //assert(m_xcoord.size()==4 and m_ycoord.size()==4);
//  m_shapefn_val.resize(4);
  m_Nis.resize(4);
  m_Nit.resize(4);

//  m_shapefn_val(0)=0.25*(1-sval)*(1-tval);
//  m_shapefn_val(1)=0.25*(1+sval)*(1-tval);
//  m_shapefn_val(2)=0.25*(1+sval)*(1+tval);
//  m_shapefn_val(3)=0.25*(1-sval)*(1+tval);

  m_Nis[0]=(tval-1)/4;
  m_Nis[1]=-(tval-1)/4;
  m_Nis[2]=(tval+1)/4;
  m_Nis[3]=-(tval+1)/4;

  m_Nit[0]= (sval-1)/4;
  m_Nit[1]=-(sval+1)/4;
  m_Nit[2]=(sval+1)/4;
  m_Nit[3]=-(sval-1)/4;

 // Eigen::Map<Eigen::VectorXd>m_x(m_xcoord.data(),m_xcoord.size());
 // Eigen::Map<Eigen::VectorXd>m_y(m_ycoord.data(),m_ycoord.size());

    Eigen::VectorXd m_x(4);
    Eigen::VectorXd m_y(4);
    for(unsigned int i=0;i<4;i++)
    {
        m_x[i]=m_xcoord[i];
        m_y[i]=m_ycoord[i];
    }

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

Eigen::MatrixXd Quad4::get_bmat()
{

    Eigen::MatrixXd bdmat(3,8);
    bdmat.setZero();

    bdmat(0,0)=m_Nis(0)*m_inverse_jacobian(0,0)+m_Nit(0)*m_inverse_jacobian(1,0);
    bdmat(0,2)=m_Nis(1)*m_inverse_jacobian(0,0)+m_Nit(1)*m_inverse_jacobian(1,0);
    bdmat(0,4)=m_Nis(2)*m_inverse_jacobian(0,0)+m_Nit(2)*m_inverse_jacobian(1,0);
    bdmat(0,6)=m_Nis(3)*m_inverse_jacobian(0,0)+m_Nit(3)*m_inverse_jacobian(1,0);

    bdmat(1,1)=m_Nis(0)*m_inverse_jacobian(0,1)+m_Nit(0)*m_inverse_jacobian(1,1);
    bdmat(1,3)=m_Nis(1)*m_inverse_jacobian(0,1)+m_Nit(1)*m_inverse_jacobian(1,1);
    bdmat(1,5)=m_Nis(2)*m_inverse_jacobian(0,1)+m_Nit(2)*m_inverse_jacobian(1,1);
    bdmat(1,7)=m_Nis(3)*m_inverse_jacobian(0,1)+m_Nit(3)*m_inverse_jacobian(1,1);

    bdmat(2,0)=m_Nis(0)*m_inverse_jacobian(0,1)+m_Nit(0)*m_inverse_jacobian(1,1);
    bdmat(2,2)=m_Nis(1)*m_inverse_jacobian(0,1)+m_Nit(1)*m_inverse_jacobian(1,1);
    bdmat(2,4)=m_Nis(2)*m_inverse_jacobian(0,1)+m_Nit(2)*m_inverse_jacobian(1,1);
    bdmat(2,6)=m_Nis(3)*m_inverse_jacobian(0,1)+m_Nit(3)*m_inverse_jacobian(1,1);

    bdmat(2,1)=m_Nis(0)*m_inverse_jacobian(0,0)+m_Nit(0)*m_inverse_jacobian(1,0);
    bdmat(2,3)=m_Nis(1)*m_inverse_jacobian(0,0)+m_Nit(1)*m_inverse_jacobian(1,0);
    bdmat(2,5)=m_Nis(2)*m_inverse_jacobian(0,0)+m_Nit(2)*m_inverse_jacobian(1,0);
    bdmat(2,7)=m_Nis(3)*m_inverse_jacobian(0,0)+m_Nit(3)*m_inverse_jacobian(1,0);

    return bdmat;


}






