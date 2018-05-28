#ifndef SHAPEFN_H
#define SHAPEFN_H
#include<memory>
#include"eigen3/Eigen/Dense"
class Shapefn
{
public:
    Shapefn(){}

    Shapefn(std::vector<double> xcoord, std::vector<double> ycoord){m_xcoord=xcoord;m_ycoord=ycoord;}
    ~Shapefn(){}
    void add_coordinates(std::vector<double> xcoord, std::vector<double> ycoord){m_xcoord=xcoord;m_ycoord=ycoord;}
    virtual void update_shapefn(double sval, double tval)=0;
    Eigen::MatrixXd get_bmat(){return m_bmat;}    //strain displacement matrix
    Eigen::MatrixXd get_jacobian(){return m_jacobian;}
    Eigen::MatrixXd& get_inverse_jacobian(){return m_inverse_jacobian;}
    Eigen::RowVectorXd get_shapefn_val(){return m_shapefn_val;}
    Eigen::RowVectorXd get_Nis(){return m_Nis;}
    Eigen::RowVectorXd get_Nit(){return m_Nit;}
    double get_detJ(){return m_detJ;}
    double get_detJ(double sval,double tval);
    double X(){}
    double Y(){}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    Eigen::MatrixXd m_bmat;
    Eigen::MatrixXd m_jacobian;
    Eigen::MatrixXd m_inverse_jacobian;
    double m_detJ;
    Eigen::RowVectorXd m_shapefn_val;
    Eigen::RowVectorXd m_Nis;
    Eigen::RowVectorXd m_Nit;
    std::vector<double> m_xcoord;
    std::vector<double> m_ycoord;

};

class Quad4:public Shapefn
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Quad4(){}
    ~Quad4(){}
    Quad4(std::vector<double> xcoord, std::vector<double> ycoord){m_xcoord=xcoord;m_ycoord=ycoord;}
    void update_shapefn( double sval,double tval);
    void set_coordinates(std::vector<double>,std::vector<double>);
    Eigen::MatrixXd get_bsmat();


};
#endif // SHAPEFN_H
