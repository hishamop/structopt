#ifndef CPS4_NLP_H
#define CPS4_NLP_H


#include "cmodel.h"
#include "coin/IpTNLP.hpp"
using namespace Ipopt;

class CPS4_NLP:public TNLP
{
public:
    CPS4_NLP(CModel* model);
    virtual ~CPS4_NLP(){}
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values){return true;}
   virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                 const IpoptData* ip_data,
                                  IpoptCalculatedQuantities* ip_cq);

    void set_nlp_info();
    void set_constraints();
    int get_Nvar(){return m_Nvar;}
    int  get_Nconstr(){return m_Nconstr;}
    int  get_Nnz_jacobian(){return m_Nnz_jacobian;}
    int  get_Nnz_hess_lagrangian();

    //dof , gradient and hessian
    void update_gradient(Eigen::RowVectorXd&);
    std::vector<double>& get_gradvals(){return m_gradvals;}

private:



  void initialize_elements_dofs();
  void update_dof(const double*);
  void update_funvals(const double*);





private:     //PRIVATE DATA
  CModel*      m_model;

  double m_objval;
  Eigen::RowVectorXd    m_fungrad;
  bool grad_set;

  int m_Nconstr;               //No of constraints
  int m_Nvar;                  //No of variables
  int m_Nnz_jacobian;          //No of nonzero in jacobian
  int m_Nnz_hess_lagrangian;   //No of nonzero in hessian of lagrgangian
//    double m_obj_val;            //obj. function value
  std::vector<double> m_xi;     //initial values
  std::vector<double> m_xi_lower_bound;
  std::vector<double> m_xi_upper_bound;
  //std::vector<double> m_g;       //Constraint equations
  std::vector<double> m_gval;  //constraint values
  std::vector<double> m_x;
  std::vector<double> m_g_val;   //Constraint value
  std::vector<double> m_gradvals; //gradient of function variables

  void set_constraint_matrix();

};



#endif // CPS4_NLP_H
