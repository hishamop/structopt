#ifndef CQUADRATRUE_H
#define CQUADRATRUE_H

#include<vector>
#include<iostream>
#include"eigen3/Eigen/Core"
class CQuadrature
{
public:
    CQuadrature(int);
    void set_quadrature(int);
    const std::vector<double> get_nodes() const
    {
        return m_nodes;
    }
    const std::vector<double> get_weight() const
    {
        return m_weight;
    }

    int get_degree() const
    {
        return m_pts;
    }

    double AreaIntegral(double(*function)(double,double));
    Eigen::MatrixXd AreaIntegral(int, double (*function)(double, double));

private:
    int                  m_pts;
    std::vector<double>  m_nodes;
    std::vector<double>  m_weight;
    void set_gauss_val(const double*,const double*);
};


#endif // CQUADRATRUE_H
