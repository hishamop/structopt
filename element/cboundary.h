#ifndef CBOUNDARY_H
#define CBOUNDARY_H

#include <vector>
#include <map>
#include<memory>
#include"eigen3/Eigen/Core"
#include"cnode.h"

using node_ptr= std::shared_ptr<CNode>;
class BoundaryFace
{
    struct TripletMatrix
    {
        std::vector<unsigned int> m_row;
        std::vector<unsigned int> m_col;
        std::vector<double> m_val;

    };

public:
    BoundaryFace();
    BoundaryFace(node_ptr n1,node_ptr n2);
    double length();    //length of face.
    double nx,ny;
    unsigned int id;
    node_ptr node1,node2;
    double pressure=0;
    int m_fixity=0;    
    int   m_index;  //index for constraint beginning.
    Eigen::MatrixXd m_constraint; //bsmat at nodes 1 and two assgined by Element class
    bool is_constraint_set;
    std::vector<double> m_boundvals;  //constraints boundary vals.
    TripletMatrix m_jac_g;  //jacobian of constraints.


};

class CBoundary
{
public:
    CBoundary(){}
    ~CBoundary(){}
    void set_constraints(const double* x,const double* g);

    std::vector<BoundaryFace> m_boundary_faces;
    std::vector<int> m_faces;
    std::vector<int> m_constrained_faces;
    std::map<unsigned int,double> m_pload;



};

class CPS4_boundary:public CBoundary
{
public:
    CPS4_boundary(){}


};

#endif // CBOUNDARY_H
