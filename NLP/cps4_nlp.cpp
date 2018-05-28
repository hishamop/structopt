
#include "cps4_nlp.h"
#include"eigen3/Eigen/Core"
#include<memory>

CPS4_NLP::CPS4_NLP(CModel *model)
    {
        m_model=model; m_Nvar=0;m_Nconstr=0; m_Nnz_jacobian=0; m_Nnz_hess_lagrangian=0;
  }

void CPS4_NLP::set_nlp_info()
{
    m_Nvar=m_model->m_nodes.size()*8 ;  //Each node has 8  DOFs.
    m_gradvals.resize(m_Nvar);
    m_Nnz_jacobian=0;

    int index=0;
    for(auto&& iter:m_model->m_boundary_elems)
    {
        m_Nconstr+=iter->m_boundary->m_boundary_faces.size()*7;
        for(auto& face:iter->m_boundary->m_boundary_faces)
        {
            face.m_index= index;
            index+=7;
        }
    }

    for(auto&&iter:m_model->m_boundary_elems)
    {

            for(auto&& it:iter->m_boundary->m_boundary_faces)
            {

                m_Nnz_jacobian+=it.m_jac_g.m_val.size();
                iter->set_NLP_constraint_matrix(it);
                m_gval.insert(m_gval.end(),it.m_boundvals.begin(),it.m_boundvals.end());
            }


    }




}

bool CPS4_NLP:: get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
      set_nlp_info();


   n=m_Nvar;
   m=m_Nconstr;
   nnz_jac_g = m_Nnz_jacobian ;
   //nnz_h_lag = 10;
   index_style = TNLP::C_STYLE;
   return true;
}


bool CPS4_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
    std::copy(m_xi_lower_bound.begin(),m_xi_lower_bound.end(),x_l);
    std::copy(m_xi_upper_bound.begin(),m_xi_upper_bound.end(),x_u);


    std::copy(m_g_val.begin(),m_g_val.end(),g_l);
    std::copy(m_g_val.begin(),m_g_val.end(),g_u);

    return true;
}

void CPS4_NLP::update_dof(const double *x)
{
    for(auto& iter:m_model->m_nodes)
    {
        iter->update_dof(x);
    }
}



/*
CPS4_NLP::set_constraints()
{
    for(auto& iter:m_model->m_boundary_elems)
    {

    }

}*/


 bool CPS4_NLP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U,
                                   Index m, bool init_lambda, Number *lambda)
 {

     for(unsigned int i=0;i<n;i++)
     {
         x[i]=0;
     }
     return true;
//         std::copy(m_xi.begin(),m_xi.end(),x);
         return true;
 }


 bool CPS4_NLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value)
 {
     if(new_x)
     {
         update_funvals(x);
     }

     obj_value = this->m_objval;
     return true;
 }



 bool CPS4_NLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f)
 {
     // All gradients are initiliased to zeros
    if(new_x)
    {
        update_funvals(x);
    }

    std::copy(m_gradvals.begin(),m_gradvals.end(),grad_f);
    return true;
 }

 bool CPS4_NLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g)
 {

     for(auto iter:m_model->m_boundary_elems)
     {
         iter->set_boundary_constraints(x,g);
     }
     return true;
 }

 bool CPS4_NLP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol, Number *values)
 {

     if (values == NULL) {
       // return the structure of the jacobian
         for(auto iter:m_model->m_boundary_elems)
         {
             int i=0;
             for(auto face:iter->m_boundary->m_boundary_faces)
             {
                 auto& row = face.m_jac_g.m_row;
                 auto& col = face.m_jac_g.m_col;
                 for(auto& rowIter:row)
                     iRow[i++] = rowIter;
                 for(auto& colIter:col)
                     jCol[i++] = colIter;

             }
         }

     }
     else {
       // return the values of the jacobian of the constraints

         for(auto& iter:m_model->m_boundary_elems)
         {
             int i=0;
             for(auto face:iter->m_boundary->m_boundary_faces)
             {
                 auto& val = face.m_jac_g.m_val;
                 for(auto& valIter:val)
                     values[i++] = valIter;
             }


         }


     }

     return true;
 }


 void CPS4_NLP::update_funvals(const double*x)
 {
     //finding objective function value and gradients
     std::fill(m_gradvals.begin(),m_gradvals.end(),0.0);
     m_objval=0;
     for(auto&iter:m_model->m_nodes)
     {
         iter->update_dof(x);
     }
     for(auto& iter:m_model->m_elements)
     {
         
         iter->calculate_fvals(m_gradvals);
         m_objval+=iter->get_objval();
     }
 }

