#ifndef FEVAL_H
#define FEVAL_H
#include "eigen3/Eigen/Core"
#include "enums.h"
#include "cquadratrue.h"



class FEval
{
public:
    //FEval(){m_quad(10);}
    FEval(Element::element_type type):m_type(type),m_quad(10),m_qpts(10){compute();}
     FEval(Element::element_type type,int quad_pts):m_type(type),m_quad(quad_pts),m_qpts(quad_pts){compute();}
    void compute();
    std::vector<Eigen::MatrixXd>& get_bsmat()
    {
        return m_bsmat;
    }
    const std::vector<double> get_qnodes(){return m_quad.get_nodes();}
    const std::vector<double> get_qweight(){return m_quad.get_weight();}
    const Eigen::MatrixXd& bsmat(unsigned int i,unsigned int j){
        return m_bsmat.at(i*m_qpts+j);
    }

//    const Eigen::MatrixXd& bdmat(std::vector<node_ptr>&,unsigned int,unsigned int);


private:
    Element::element_type m_type;
    CQuadrature  m_quad;
    int m_qpts;
    std::vector<Eigen::MatrixXd> m_bsmat;
//    std::vector<Eigen::MatrixXd> m_bdmat;

    void cps4_bsmat();
    void cps4_bdmat();
};

#endif // FEVAL_H
