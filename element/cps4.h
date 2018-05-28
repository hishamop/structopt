#ifndef CPS4_H
#define CPS4_H
#include "celement.h"
#include "cnode.h"
#include <vector>
#include<memory>
#include"eigen3/Eigen/Core"
#include"material/cmaterial.h"
#include"feval.h"

class CPS4:public CElement
{
public:
    CPS4(){m_id=0; m_type=Element::CPS4; m_nodes.reserve(4);}
    CPS4(int id);
    CPS4(int id,std::vector<node_ptr> nlist);
    ~CPS4(){}

    //virtual void get_fval   (const double* x,CMaterial*);
    //virtual void set_grad_f (const double* x,const double* g, CMaterial*){}
    //void set_node_dof(){}

    virtual     void set_boundary_info();
    void set_boundary_faces();
    virtual void set_constrained_faces();




    virtual void set_boundary_constraints(const double*, double*);
   // virtual void evaluate_constraints();
//    virtual void set_jacobian_structure(const double*,const double*);

    std::vector<node_ptr> get_face(unsigned int) const ;
    int elements_shared_by_face    (unsigned int ) const;
    Element::element_type get_element_type() const;
    int faces_count() const;
    bool is_boundary_element() const;
    unsigned int getindex() const;




    //transformation and other matrices
    void set_local_stress_dof(Eigen::VectorXd&);
    Eigen::MatrixXd get_dmat();
    //Eigen::Vector3d stress_local(double, double);
    Eigen::Vector3d strain_local(double,double);


    std::vector<double> get_stress_vector(double sval, double tval, const double* dof) const;
   // double get_objval(const std::vector<double>& x, material_ptr material, FEval& vals) const;
    std::vector<double> get_ref_dof(std::vector<double>& Dof) const;
    std::vector<double> get_standard_stress_dof(const double* var_list);
    std::vector<double> get_standard_stress_dof(const std::vector<double>& var_list);
    Eigen::VectorXd     get_standard_stress_dof(Eigen::VectorXd& var);
    void calculate_fvals(std::vector<double>&);
    void set_NLP_constraint_matrix(BoundaryFace&);
    void set_constraints_boundvals();
 //   std::vector<unsigned int> get_node_ids();

    //shape functions
    Eigen::Matrix2d jacobian_matrix(double sval, double tval);
    Eigen::Matrix2d inverse_jacobian(double sval,double tval);
    Eigen::Matrix<double,3,8> get_bmat(double sval,double tval);
    Eigen::MatrixXd get_bsmat_local(double x,double y);
    double detJ(double,double);

    int get_dof_count() const{return 32;}

    void add_shape_function();


private:

    static FEval m_val;
    static CMaterial* m_material;
};


#endif // CPS4_H
