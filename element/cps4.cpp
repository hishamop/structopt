#include "cps4.h"
#include"material/cmaterial.h"
#include"cquadratrue.h"
#include"feval.h"
#include <cassert>
#include <algorithm>
#include <iterator>
#include<cmath>
#include<memory>
#include"shape/shapefn.h"
FEval CPS4::m_val(Element::CPS4,10);
  //static object initializer
CMaterial* CPS4::m_material=new STEEL();

CPS4::CPS4(int id)
{m_id=id; m_type=Element::CPS4;m_nodes.reserve(4);
    m_shape=std::make_shared<Quad4>();

}

void CPS4::add_shape_function()
{
    m_shape=std::make_shared<Quad4>(this->get_xcoords(),this->get_ycoords());

}

CPS4::CPS4(int id, std::vector<node_ptr> nodes)
{
    m_type =Element::CPS4;
    m_id=id;
    m_nodes.reserve(4);
    m_nodes=nodes;
    m_shape=std::make_shared<Quad4>(this->get_xcoords(),this->get_ycoords());

}





std::vector<node_ptr> CPS4::get_face(unsigned int i) const
{
    std::vector<node_ptr> temp;
    if(i !=3)
    {
        temp.push_back(m_nodes.at(i));
        temp.push_back(m_nodes.at(i+1));
    }
    else
    {
        temp.push_back(m_nodes.at(3));
        temp.push_back(m_nodes.at(0));
    }

    return temp;
}

int CPS4::elements_shared_by_face(unsigned int i) const
{
    std::vector<node_ptr> face = get_face(i);
    std::vector<unsigned int> v1 = face.at(0)->get_shared_elements_id() ; std::sort(v1.begin(),v1.end());
    std::vector<unsigned int> v2 = face.at(1)->get_shared_elements_id() ; std::sort(v2.begin(),v2.end());
    std::vector<unsigned int> intersection ;
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::inserter(intersection,intersection.begin()));
    return intersection.size();
}

Element::element_type CPS4::get_element_type() const
{
    return m_type;
}



int CPS4::faces_count() const
{
    return 4;
}

bool CPS4::is_boundary_element( ) const
{
    bool flag =false;
    for(int i=0; i<4; i++)
    {
        if(elements_shared_by_face(i) == 1)
        {
            flag = true;
            break;
        }
    }
    return flag;
}

unsigned int CPS4::getindex() const
{
    return m_id;
}


/*std::vector<double>    cps4::get_stress_vector(double sval, double tval, std::vector<double> dof,
                                                     cell_type = m_type
                                                     )                   const
{


    std::vector<double> stress_vec(3);
    std::vector<double> sxx_coeff(24,0.0);
    std::vector<double> syy_coeff(24,0.0);
    std::vector<double> sxy_coeff(24,0.0);

    sxx_coeff[0]=  (3*tval*(sval - 1)*(sval^2 + sval + 5*tval^2 - 5))/8;
    sxx_coeff[1]= 3*tval*(sval - 1)^2*(sval + 1)/8;
    sxx_coeff[2] = (sval - 1)*(3*sval^2*tval - sval^2 + 3*sval*tval - sval + 15*tval^3 - 3*tval^2 - 15*tval + 3)/8;
    sxx_coeff[4] = (sval - 1)*(tval - 1)*(5*tval^2 + 2*tval - 1)/8;
    sxx_coeff[5] = (3*tval - 1)*(sval - 1)^2*(sval + 1)/8;
    sxx_coeff[6]= 3*tval*(sval + 1)*(- sval^2 + sval - 5*tval^2 + 5)/8;
    sxx_coeff[7] = 3*tval*(sval - 1)*(sval + 1)^2/8;
    sxx_coeff[8] = -(sval + 1)*(3*sval^2*tval - sval^2 - 3*sval*tval + sval + 15*tval^3 - 3*tval^2 - 15*tval + 3)/8;
    sxx_coeff[10]= -(sval + 1)*(tval - 1)*(5*tval^2 + 2*tval - 1)/8;
    sxx_coeff[11] = (3*tval - 1)*(sval - 1)*(sval + 1)^2/8;
    sxx_coeff[12] = 3*tval*(sval + 1)*(sval^2 - sval + 5*tval^2 - 5)/8;
    sxx_coeff[13] = -3*tval*(sval - 1)*(sval + 1)^2/8;
    sxx_coeff[14] = (sval + 1)*(-3*sval^2*tval - sval^2 + 3*sval*tval + sval - 15*tval^3 - 3*tval^2 + 15*tval + 3)/8;
    sxx_coeff[16] = (sval + 1)*(tval + 1)*(5*tval^2 - 2*tval - 1)/8;
    sxx_coeff[17] = (3*tval + 1)*(sval - 1)*(sval + 1)^2/8;
    sxx_coeff[18] = -(3*tval*(sval - 1)*(sval^2 + sval + 5*tval^2 - 5))/8;
    sxx_coeff[19] = -3*tval*(sval - 1)^2*(sval + 1)/8;
    sxx_coeff[20]= (sval - 1)*(3*sval^2*tval + sval^2 + 3*sval*tval + sval + 15*tval^3 + 3*tval^2 - 15*tval - 3)/8;
    sxx_coeff[22] = (sval - 1)*(tval + 1)*(- 5*tval^2 + 2*tval + 1)/8;
    sxx_coeff[23] = (3*tval + 1)*(sval - 1)^2*(sval + 1)/8;


    syy_coeff[0] = (3*sval*(tval - 1)*(5*sval^2 + tval^2 + tval - 5))/8;
    syy_coeff[1] = (tval - 1)*(15*sval^3 - 3*sval^2 + 3*sval*tval^2 + 3*sval*tval - 15*sval - tval^2 - tval + 3)/8;
    syy_coeff[2] = 3*sval*(tval - 1)^2*(tval + 1)/8;
    syy_coeff[3] = (sval - 1)*(tval - 1)*(5*sval^2 + 2*sval - 1)/8;
    syy_coeff[5] = (3*sval - 1)*(tval - 1)^2*(tval + 1)/8;
    syy_coeff[6] = -3*sval*(tval - 1)*(5*sval^2 + tval^2 + tval - 5)/8;
    syy_coeff[7]= (tval - 1)*(15*sval^3 + 3*sval^2 + 3*sval*tval^2 + 3*sval*tval - 15*sval + tval^2 + tval - 3)/8;
    syy_coeff[8]= -3*sval*(tval - 1)^2*(tval + 1)/8;
    syy_coeff[9]= ((sval + 1)*(tval - 1)*(- 5*sval^2 + 2*sval + 1))/8;
    syy_coeff[11] = (3*sval + 1)*(tval - 1)^2*(tval + 1)/8;
    syy_coeff[12] = 3*sval*(tval + 1)*(5*sval^2 + tval^2 - tval - 5)/8;
    syy_coeff[13S] = (tval + 1)*(-15*sval^3 - 3*sval^2 - 3*sval*tval^2 + 3*sval*tval + 15*sval - tval^2 + tval + 3)/8;
    syy_coeff[14]= -3*sval*(tval - 1)*(tval + 1)^2/8;
    syy_coeff[15] = (sval + 1)*(tval + 1)*(5*sval^2 - 2*sval - 1)/8;
    syy_coeff[17]= (3*sval + 1)*(tval - 1)*(tval + 1)^2/8;
    syy_coeff[18] = (3*sval*(tval + 1)*(-5*sval^2 - tval^2 + tval + 5))/8;
    syy_coeff[19] = (tval + 1)*(-15*sval^3 + 3*sval^2 - 3*sval*tval^2 + 3*sval*tval + 15*sval + tval^2 - tval - 3)/8;
    syy_coeff[20] = 3*sval*(tval - 1)*(tval + 1)^2/8;
    syy_coeff[21] = -(sval - 1)*(tval + 1)*(5*sval^2 + 2*sval - 1)/8;
    syy_coeff[23] = (3*sval - 1)*(tval - 1)*(tval + 1)^2/8;

    sxy_coeff[0] = -(15*sval^4 + 18*sval^2*tval^2 - 36*sval^2 + 15*tval^4 - 36*tval^2 + 24)/32;
    sxy_coeff[1]= -(3*sval + 1)*(sval - 1)*(5*sval^2 + 2*sval + 6*tval^2 - 9)/32;
    sxy_coeff[2] = -(3*tval + 1)*(tval - 1)*(6*sval^2 + 5*tval^2 + 2*tval - 9)/32;
    sxy_coeff[3] = -(5*sval + 1)*(sval - 1)^2*(sval + 1)/32;
    sxy_coeff[4] = -(5*tval + 1)*(tval - 1)^2*(tval + 1)/32;
    sxy_coeff[5] = -(3*sval + 1)*(3*tval + 1)*(sval - 1)*(tval - 1)/16;
    sxy_coeff[6] = (15*sval^4 + 18*sval^2*tval^2 - 36*sval^2 + 15*tval^4 - 36*tval^2 + 24)/32;
    sxy_coeff[7] = (3*sval - 1)*(sval + 1)*(- 5*sval^2 + 2*sval - 6*tval^2 + 9)/32;
    sxy_coeff[8] = (3*tval + 1)*(tval - 1)*(6*sval^2 + 5*tval^2 + 2*tval - 9)/32;
    sxy_coeff[9] = (5*sval - 1)*(sval - 1)*(sval + 1)^2/32;
    sxy_coeff[10]= (5*tval + 1)*(tval - 1)^2*(tval + 1)/32;
    sxy_coeff[11] = -(3*sval - 1)*(3*tval + 1)*(sval + 1)*(tval - 1)/16;
    sxy_coeff[12]= -(15*sval^4 + 18*sval^2*tval^2 - 36*sval^2 + 15*tval^4 - 36*tval^2 + 24)/32;
    sxy_coeff[13] = (3*sval - 1)*(sval + 1)*(5*sval^2 - 2*sval + 6*tval^2 - 9)/32;
    sxy_coeff[14]= (3*tval - 1)*(tval + 1)*(6*sval^2 + 5*tval^2 - 2*tval - 9)/32;
    sxy_coeff[15] = -(5*sval - 1)*(sval - 1)*(sval + 1)^2/32;
    sxy_coeff[16]= -(5*tval - 1)*(tval - 1)*(tval + 1)^2/32;
    sxy_coeff[17] = -(3*sval - 1)*(3*tval - 1)*(sval + 1)*(tval + 1)/16;
    sxy_coeff[18] = (15*sval^4 + 18*sval^2*tval^2 - 36*sval^2 + 15*tval^4 - 36*tval^2 + 24)/32;
    sxy_coeff[19] = (3*sval + 1)*(sval - 1)*(5*sval^2 + 2*sval + 6*tval^2 - 9)/32;
    sxy_coeff[20] = (3*tval - 1)*(tval + 1)*(-6*sval^2 - 5*tval^2 + 2*tval + 9)/32;
    sxy_coeff[21] = (5*sval + 1)*(sval - 1)^2*(sval + 1)/32;
    sxy_coeff[22] = (5*tval - 1)*(tval - 1)*(tval + 1)^2/32;
    sxy_coeff[23] = -(3*sval + 1)*(3*tval - 1)*(sval - 1)*(tval + 1)/16;


    return stress_vec;

}
*/
/*double CPS4::get_objval(const std::vector<double>& x, material_ptr material, FEval &val) const
{

    double objval;
    auto bsmat = val.get_bsmat();
    //using the cached bsmat matrices.
    for (auto iter:bsmat)
    {

    }

    return objval;

}*/

void CPS4::set_local_stress_dof(Eigen::VectorXd& dof)
{

    double dxds,dyds,dxdt,dydt;
    Eigen::MatrixXd jac=m_shape->get_jacobian();
    dxds=jac(0,0);
    dxdt=jac(0,1);
    dyds=jac(1,0);
    dydt=jac(1,1);

    int i=0;
    for(auto&iter:m_nodes)
    {
       const double* sdof=iter->get_stress_dof();
       dof[i*6]  = sdof[0];
       dof[i*6+1]= dxds*sdof[1]+dyds*sdof[2];
       dof[i*6+2]= dxdt*sdof[1]+dydt*sdof[2];
       dof[i*6+3]= dxds*dxds*sdof[3]+dyds*dyds*sdof[4]+2*dxds*dyds*sdof[5];
       dof[i*6+4]= dxdt*dxdt*sdof[3]+dydt*dydt*sdof[4]+2*dxdt*dydt*sdof[5];
       dof[i*6+5]= dxds*dxdt*sdof[3]+dyds*dydt*sdof[4]+(dxds*dydt+dxdt*dyds)*sdof[5];;
       i++;
    }
}


void CPS4::set_boundary_info()
{
    for(int faceId=0;faceId<4;faceId++)
    {
        std::vector<node_ptr> face = get_face(faceId);
        std::vector<unsigned int> v1 = face.at(0)->get_shared_elements_id() ; std::sort(v1.begin(),v1.end());
        std::vector<unsigned int> v2 = face.at(1)->get_shared_elements_id() ; std::sort(v2.begin(),v2.end());
        std::vector<unsigned int> intersection ;
        set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::inserter(intersection,intersection.begin()));

        if(intersection.size() == 1)
        {
            is_boundary = true;
            face.at(0)->on_boundary();
            face.at(1)->on_boundary();
            BoundaryFace temp(face[0],face[1]);
            temp.id=faceId;
            m_boundary->m_boundary_faces.push_back(temp);
            m_boundary->m_faces.push_back(faceId);
        }
    }
}

void CPS4::set_boundary_faces()
{
    for( int i=0; i<4; i++)
    {
        if(elements_shared_by_face(i) == 1)
        {
            m_boundary->m_faces.push_back(i);
            m_boundary->m_pload[i]=0.00;
//            assign_traction(i);
        }

    }
}

void CPS4::set_constrained_faces()
{
    if(m_nodes[0]->is_constrined() && m_nodes[1]->is_constrined())
    {
        m_boundary->m_constrained_faces.push_back(0);
        m_boundary->m_pload.erase(0);
    }
    if(m_nodes[1]->is_constrined() && m_nodes[2]->is_constrined())
    {
        m_boundary->m_constrained_faces.push_back(1);
        m_boundary->m_pload.erase(1);
    }

    if(m_nodes[2]->is_constrined() && m_nodes[3]->is_constrined())
    {
        m_boundary->m_constrained_faces.push_back(2);
        m_boundary->m_pload.erase(2);
    }
    if(m_nodes[3]->is_constrined() && m_nodes[0]->is_constrined())
    {
        m_boundary->m_constrained_faces.push_back(3);
        m_boundary->m_pload.erase(3);
    }

}



void CPS4::calculate_fvals(std::vector<double>& gval)
{
    const std::vector<double>& qn = m_val.get_qnodes();
    const std::vector<double>& qw = m_val.get_qweight();

    double objval =0.0;
  //  Eigen::RowVectorXd sdf_grad(24);
    Eigen::RowVectorXd disp_grad = Eigen::VectorXd::Zero(8);
    Eigen::VectorXd sdf_grad= Eigen::VectorXd::Zero(24);
    Eigen::VectorXd stress_dof = get_stress_dof();
    Eigen::VectorXd strain_dof = get_strain_dof();


    for(unsigned int i=0;i!=qn.size(); i++)
    {
        for(unsigned int j=0; j!=qn.size();j++)
        {
            double sval = qn[i];
            double tval = qn[j];
            if(!m_shape->m_set_coordinate)
                 m_shape->add_coordinates(this->get_xcoords(),this->get_ycoords());
            m_shape->update_shapefn(sval,tval);
            Eigen::MatrixXd Bs = get_bsmat_local(i,j);
            Eigen::Vector3d stress = Bs*stress_dof;
            Eigen::MatrixXd Bd = get_bmat(sval,tval);
            Eigen::Vector3d strain = Bd* strain_dof;

//            Eigen::Vector3d strain = strain_local(sval,tval);
            Eigen::Vector3d cnfn = m_material->get_constitutive_function(stress,strain);
            auto WxJ  = qw[i] * qw[j] * m_shape->get_detJ();
//            double delval= cnfn.squaredNorm();
            double delval=cnfn.norm();
            objval+=delval*WxJ;

            //SETTING GRADIENT OF STRESS DOF

            Eigen::MatrixXd CxBs = m_material->get_cmat()*Bs;
//            sdf_grad+= -2*CxBs.transpose()*cnfn*WxJ;
            Eigen::MatrixXd Bstrans = CxBs.transpose();
            Eigen::VectorXd gradval= (-WxJ/delval)*Bstrans*cnfn;
            sdf_grad=sdf_grad+ gradval ;
//            sdf_grad+= -1*CxBs.transpose()*cnfn*WxJ/delval;


            //SETTING GRADIENT OF DISP.DOF
            Eigen::MatrixXd bmat=m_shape->get_bmat();
            disp_grad+=(WxJ/delval)*bmat.transpose()*cnfn;
//            disp_grad+= bmat.transpose()*cnfn*2.0*WxJ;;
        }
   }


    auto node_ids=get_node_ids(); unsigned int i=0;
    for(auto&iter:node_ids)
    {

        unsigned int pos=(iter-1)*8;
        gval[pos]+=sdf_grad[i*6];
        gval[pos+1]+=sdf_grad[i*6+1];
        gval[pos+2]+=sdf_grad[i*6+2];
        gval[pos+3]+=sdf_grad[i*6+3];
        gval[pos+4]+=sdf_grad[i*6+4];
        gval[pos+5]+=sdf_grad[i*6+5];

        gval[pos+6]+=disp_grad[i*2];
        gval[pos+7]+=disp_grad[i*2+1];

        i++;
    }

    m_objval=objval;
}

Eigen::MatrixXd CPS4::get_dmat()
{
    Eigen::MatrixXd dmat(3,3);
    double dtdy,dsdy,dtdx,dsdx;
    Eigen::MatrixXd& invjac= m_shape->get_inverse_jacobian();
    dsdx= invjac(0,0);
    dsdy= invjac(0,1);
    dtdx= invjac(1,0);
    dtdy=invjac(1,1);

    dmat(0,0)=dtdy*dtdy;
    dmat(0,1)=dsdy*dsdy;
    dmat(0,2)=-2*dsdy*dtdy;
    dmat(1,0)=dtdx*dtdx;
    dmat(1,1)=dsdx*dsdx;
    dmat(1,2)=-2*dsdx*dtdx;
    dmat(2,0)=-dtdx*dtdy;
    dmat(2,1)=-dsdx*dsdy;
    dmat(2,2)=dsdx*dtdy+dsdy*dtdx;
    return dmat;
}

Eigen::Vector3d CPS4::strain_local(double sval, double tval)
{

    auto bdmat=get_bmat(sval,tval);
    Eigen::VectorXd uvec=get_strain_dof();
    return bdmat*uvec;
}

Eigen::Matrix<double,3,8> CPS4::get_bmat(double,double)
{
    Eigen::MatrixXd inv= m_shape->get_inverse_jacobian();
    Eigen::VectorXd Nis,Nit;
    Nis=m_shape->get_Nis();
    Nit=m_shape->get_Nit();
    Eigen::MatrixXd bdmat(3,8);
    bdmat.setZero();

    bdmat(0,0)=Nis(0)*inv(0,0)+Nit(0)*inv(1,0);
    bdmat(0,2)=Nis(1)*inv(0,0)+Nit(1)*inv(1,0);
    bdmat(0,4)=Nis(2)*inv(0,0)+Nit(2)*inv(1,0);
    bdmat(0,6)=Nis(3)*inv(0,0)+Nit(3)*inv(1,0);

    bdmat(1,1)=Nis(0)*inv(0,1)+Nit(0)*inv(1,1);
    bdmat(1,3)=Nis(1)*inv(0,1)+Nit(1)*inv(1,1);
    bdmat(1,5)=Nis(2)*inv(0,1)+Nit(2)*inv(1,1);
    bdmat(1,7)=Nis(3)*inv(0,1)+Nit(3)*inv(1,1);

    bdmat(2,0)=Nis(0)*inv(0,1)+Nit(0)*inv(1,1);
    bdmat(2,2)=Nis(1)*inv(0,1)+Nit(1)*inv(1,1);
    bdmat(2,4)=Nis(2)*inv(0,1)+Nit(2)*inv(1,1);
    bdmat(2,6)=Nis(3)*inv(0,1)+Nit(3)*inv(1,1);

    bdmat(2,1)=Nis(0)*inv(0,0)+Nit(0)*inv(1,0);
    bdmat(2,3)=Nis(1)*inv(0,0)+Nit(1)*inv(1,0);
    bdmat(2,5)=Nis(2)*inv(0,0)+Nit(2)*inv(1,0);
    bdmat(2,7)=Nis(3)*inv(0,0)+Nit(3)*inv(1,0);

    return bdmat;

}

void CPS4:: set_boundary_constraints(const double* x, double* g)
{
    Eigen::VectorXd phi = get_stress_dof();
    for(auto& face:m_boundary->m_boundary_faces)
    {
        assert(face.m_index!= -1); // check whether indices assigned
        auto i= face.m_index;
        g[i]   = face.node1->sdof(4)*face.nx - face.node1->sdof(5)*face.ny;
        g[i+1] = face.node1->sdof(3)*face.ny - face.node1->sdof(5)*face.nx;
        g[i+2] = face.node2->sdof(4)*face.nx - face.node2->sdof(5)*face.ny;
        g[i+3] = face.node2->sdof(3)*face.ny - face.node2->sdof(5)*face.nx;

        g[i+4] = face.m_constraint.row(0)*phi;
        g[i+5] = face.m_constraint.row(1)*phi;
        g[i+6] = face.m_constraint.row(2)*phi;

        if(face.m_fixity !=0)
        {
            auto df1 = face.node1->get_displacement_dof();
            auto df2 = face.node2->get_displacement_dof();
            g[i]=g[i]-df1[0];
            g[i+1]=g[i+1]-df1[1];
            g[i+2]=g[i+2]-df2[0];
            g[i+3]=g[i+3]-df2[1];

            //force boundary condition
            g[i+4]=g[i+4]- df1[0]*face.length();
            g[i+5]=g[i+5]- df1[1]*face.length();

            //moment boundary condition
            double x1= face.node1->x();
            double x2= face.node2->x();
            double y1= face.node1->y();
            double y2= face.node2->y();
            g[i+6]=g[i+6]+(x2*x2-x1*x1)*df1[1]/2 - (y2*y2-y1*y1)*df1[0]/2;
        }
    }
}

void CPS4::set_NLP_constraint_matrix()
{
for(auto& face:m_boundary->m_boundary_faces)
 {
    if(!m_shape->m_set_coordinate)
    {
        auto xcoord=this->get_xcoords();
        auto ycoords=this->get_ycoords();
        m_shape->add_coordinates(xcoord,ycoords);
    }

    auto qn=m_val.get_qnodes();
    auto qw=m_val.get_qweight();

    //jacobian of constraints
    if(face.m_fixity==0)
    {
        auto n1=face.node1->get_id();
        auto n2=face.node2->get_id();

        auto& row=face.m_jac_g.m_row;
        auto& col= face.m_jac_g.m_col;
        auto& val=face.m_jac_g.m_val;       
        row.push_back(0);   col.push_back((n1-1)*8+4); val.push_back(face.nx);
        row.push_back(0);   col.push_back((n1-1)*8+5); val.push_back(-face.ny);
        row.push_back(1);   col.push_back((n1-1)*8+3); val.push_back(face.ny);
        row.push_back(1);   col.push_back((n1-1)*8+5); val.push_back(-face.nx);
        row.push_back(2);   col.push_back((n2-1)*8+4); val.push_back(face.nx);
        row.push_back(2);   col.push_back((n2-1)*8+5); val.push_back(-face.ny);
        row.push_back(3);   col.push_back((n2-1)*8+3); val.push_back(face.ny);
        row.push_back(3);   col.push_back((n2-1)*8+5); val.push_back(-face.nx);

    }

    if(face.m_fixity!=0)
    {
        auto n1=face.node1->get_id();
        auto n2=face.node1->get_id();

        auto& row=face.m_jac_g.m_row;
        auto& col= face.m_jac_g.m_col;
        auto& val=face.m_jac_g.m_val;

        row.push_back(0);   col.push_back((n1-1)*8+4); val.push_back(face.nx);
        row.push_back(0);   col.push_back((n1-1)*8+5); val.push_back(-face.ny);
        row.push_back(0);   col.push_back((n1-1)*8+6); val.push_back(-1);


        row.push_back(1);   col.push_back((n1-1)*8+3); val.push_back(face.ny);
        row.push_back(1);   col.push_back((n1-1)*8+5); val.push_back(-face.nx);
        row.push_back(1);   col.push_back((n1-1)*8+7); val.push_back(-1);


        row.push_back(2);   col.push_back((n1-1)*8+4); val.push_back(face.nx);
        row.push_back(2);   col.push_back((n1-1)*8+5); val.push_back(-face.ny);
        row.push_back(2);   col.push_back((n1-1)*8+6); val.push_back(-1);

        row.push_back(3);   col.push_back((n1-1)*8+3); val.push_back(face.ny);
        row.push_back(3);   col.push_back((n1-1)*8+5); val.push_back(-face.nx);
        row.push_back(3);   col.push_back((n1-1)*8+7); val.push_back(-1);


    }


    Eigen::Matrix<double,3,24> Kmat;Kmat.setZero();






    Eigen::Matrix<double,1,24> Fx,Fy,M;
    Fx.setZero();Fy.setZero();M.setZero();
    for(unsigned int i=0;i<qn.size();i++)
    {
        for(unsigned int j=0;j<qn.size();j++)
        {


            double sval=qn[i];
            double tval=qn[j];

            //m_shape->add_coordinates(this->get_xcoords(),this->get_ycoords());

            m_shape->update_shapefn(sval,tval);

            Eigen::MatrixXd J = m_shape->get_jacobian();

            double dxds,dxdt,dyds,dydt;
            dxds=J(0,0);dxdt=J(0,1);dyds=J(1,0);dydt=J(1,1);
            double Jds= sqrt(dxds*dxds+dyds*dyds);
            double Jdt= sqrt(dxdt*dxdt+dydt*dydt);

            Eigen::MatrixXd Bs =get_bsmat_local(i,j);
             Eigen::Matrix<double,1,24> Sx,Sy,Sxy,Tx,Ty;
             Sx.setZero();Sy.setZero();Sxy.setZero();Tx.setZero();Ty.setZero();

             Sx=Bs.row(1);
             Sy=Bs.row(0);
             Sxy=Bs.row(2);

             Tx=face.nx*Sx+face.ny*Sxy;
             Ty=face.ny*Sy+face.nx*Sxy;



             if(face.id==0 )
             {
                 Fx.noalias() += qw[i]*qw[j]*Tx*Jds;
                 Fy.noalias() += qw[i]*qw[j]*Ty*Jds;
                 M.noalias()  +=  qw[i]*qw[j]*(m_shape->X()*Ty-m_shape->Y()*Tx)*Jds;
             }

             if(face.id==1)
             {
                 Fx.noalias() += qw[i]*qw[j]*Tx*Jdt;
                 Fy.noalias() += qw[i]*qw[j]*Ty*Jdt;
                 M.noalias()  +=  qw[i]*qw[j]*(m_shape->X()*Ty-m_shape->Y()*Tx)*Jdt;
             }

            if(face.id==2)
             {
                 Fx.noalias() += 1*qw[i]*qw[j]*Tx*Jds;
                 Fy.noalias() += 1*qw[i]*qw[j]*Ty*Jds;
                 M.noalias()  +=  1*qw[i]*qw[j]*(m_shape->X()*Ty-m_shape->Y()*Tx)*Jds;
             }

             if(face.id==3)
             {
                 Fx.noalias() += -1*qw[i]*qw[j]*Tx*Jdt;
                 Fy.noalias() += -1*qw[i]*qw[j]*Ty*Jdt;
                 M.noalias()  += -1*qw[i]*qw[j]*(m_shape->X()*Ty-m_shape->Y()*Tx)*Jdt;
             }


        }
    }

    Kmat.row(0) = Fx;
    Kmat.row(1) = Fy;
    Kmat.row(2)  =M;

    face.m_constraint = Kmat;

    if(face.m_fixity==0)
    {
         unsigned int valN= -1; // increases from 0to23 inthe loop.
        std::vector<unsigned int> n=this->get_node_ids();
        for(auto id:n)
        {
            auto& row=face.m_jac_g.m_row;
            auto& col= face.m_jac_g.m_col;
            auto& val=face.m_jac_g.m_val;

            for(unsigned int colIndex= (id-1)*8; colIndex<(id-1)*8+6;colIndex++)
            {
               valN++;
                row.push_back(4); col.push_back(colIndex); val.push_back(Fx[valN]);
                row.push_back(5); col.push_back(colIndex); val.push_back(Fy[valN]);
                row.push_back(6); col.push_back(colIndex); val.push_back(M[valN]);
            }

        }
    }

    if(face.m_fixity!=0)
    {

        auto n=this->get_node_ids();
        auto& row=face.m_jac_g.m_row;
        auto& col= face.m_jac_g.m_col;
        auto& val=face.m_jac_g.m_val;
        auto Tx_index=(face.node1->get_id()-1)*8+6;
        auto Ty_index=(face.node1->get_id()-1)*8+7;
        for(auto id:n)
        {
            unsigned int valN=-1;
            for(unsigned int colIndex= (id-1)*8; colIndex<(id-1)*8+6;colIndex++)
            {
                valN++;
                row.push_back(4);  col.push_back(colIndex); val.push_back(Fx[valN]);
                row.push_back(5);  col.push_back(colIndex); val.push_back(Fy[valN]);
                row.push_back(6);  col.push_back(colIndex); val.push_back(M[valN]);
            }
        }

        row.push_back(4); col.push_back(Tx_index); val.push_back(-face.length());
        row.push_back(5); col.push_back(Ty_index); val.push_back(-face.length());

        double x1= face.node1->x();
        double x2= face.node2->x();
        double y1= face.node1->y();
        double y2= face.node2->y();
        //g[i+6]=g[i+6]-(x2*x2-x1*x1)*df1[1]/2 - (y2*y2-y1*y1)*df1[0]/2;
        row.push_back(6); col.push_back(Ty_index); val.push_back((x2*x2-x1*x1)/2);
        row.push_back(6); col.push_back(Tx_index); val.push_back(-(y2*y2-y1*y1)/2);

    }

    face.is_constraint_set=true;
 }
}




//Set by CPS4 Element
void CPS4::set_constraints_boundvals()
{
    std::vector<double> m_bval;
    m_bval.resize(7);
    for(auto& face:m_boundary->m_boundary_faces)
    {
        if(face.pressure==0) {face.m_boundvals=m_bval;continue;}
        if(face.m_fixity==0)
        {
            //direction cosines of normals.
            double Nx = face.ny; double Ny=-face.nx;
            double Tx,Ty;
           Tx=  m_bval[0] =m_bval[2]= face.pressure*Nx;;  //Tx
           Ty = m_bval[1]= m_bval[3]= face.pressure*Ny;   //Ty

            m_bval[4] = Tx*face.length();
            m_bval[5] = Ty*face.length();

            //integral of xTy - yTx  --> x^2*Ty/2-y^2*Tx/2
            //limit of integrals
                    double x1= face.node1->x();
                    double x2= face.node2->x();
                    double y1= face.node1->y();
                    double y2= face.node2->y();

            m_bval[6] = (x2*x2-x1*x1)*Ty/2 - (y2*y2-y1*y1)*Tx/2;
        }
        if(face.m_fixity!=0)
        {
            m_bval.resize(7); // set to zero
        }

    face.m_boundvals=m_bval;

    }
}

Eigen::VectorXd CPS4::get_stress_dof()
{
    Eigen::VectorXd phi(24);
    for(unsigned int i=0;i<4;i++)
    {
       auto sdf=m_nodes[i]->get_stress_dof();
       for(unsigned int j=0;j<6;j++)
           phi[i*6+j]=sdf[j];
    }
    return phi;
}

/*
std::vector<unsigned int>CPS4::get_node_ids()
{
  /*  std::vector<unsigned int> Ids;
    Ids.clear();
    Ids.reserve(4);
    auto x=m_nodes[0]->get_id();
    Ids.push_back(x);
    Ids.push_back(m_nodes[1]->get_id());
    Ids.push_back(m_nodes[2]->get_id());
    Ids.push_back(m_nodes[3]->get_id());



    std::vector<unsigned int> nodeIds(4);

    nodeIds[0]=m_nodes[0]->get_id();
    nodeIds[1]=m_nodes[1]->get_id();
    nodeIds[2]=m_nodes[2]->get_id();
    nodeIds[3]=m_nodes[3]->get_id();
    return nodeIds;
}*/

Eigen::MatrixXd CPS4::get_bsmat_local(unsigned int i,unsigned int j)
{
    Eigen::MatrixXd J = m_shape->get_jacobian();
    Eigen::MatrixXd D= get_dmat();
    double dxds,dxdt,dyds,dydt;
   dxds=J(0,0);dxdt=J(0,1);dyds=J(1,0);dydt=J(1,1);
    Eigen::Matrix<double,6,6> tmat;
    tmat<<  1,0,0,0,0,0,
            0,dxds,dyds,0,0,0,
            0,dxdt,dydt,0,0,0,
            0,0,0,dxds*dxds,dyds*dyds,2*dxds*dyds,
            0,0,0,dxdt*dxdt,dydt*dydt,2*dxdt*dydt,
            0,0,0,dxds*dxdt,dyds*dydt,dxds*dydt+dxdt*dyds;
    Eigen::Matrix<double,6,6>zero;zero.setZero();

    Eigen::Matrix<double,24,24> Tmat;
    Tmat<<tmat,zero,zero,zero,
            zero,tmat,zero,zero,
            zero,zero,tmat,zero,

            zero,zero,zero,tmat;
    Eigen::MatrixXd Bs = D*m_val.bsmat(i,j)*Tmat;
    return Bs;

}

Eigen::VectorXd CPS4::get_strain_dof()
{
    Eigen::VectorXd uvec(8);
    int i=0;
    for(auto&iter:m_nodes)
    {
        auto u= iter->get_displacement_dof();
        uvec(2*i)=u[0];
        uvec(2*i+1)=u[1];
        i++;
    }

    return uvec;
}
