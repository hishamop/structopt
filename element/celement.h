#ifndef CELEMENT_H
#define CELEMENT_H
#include "material/cmaterial.h"
//#include "shape/shapefn.h"
#include "cnode.h"
#include"enums.h"
#include "feval.h"
#include "cboundary.h"
#include <vector>
#include <memory>
#include"shape/shapefn.h"
using boundary_ptr= std::shared_ptr<CBoundary>;
using shape_ptr = std::shared_ptr<Shapefn>;
using node_ptr=std::shared_ptr<CNode>;
using material_ptr= std::shared_ptr<CMaterial>;
class FEVal;
class CElement
{


public:
    CElement(){m_shape=nullptr;}
    CElement(int id){m_id=id; m_shape=nullptr;}
    CElement(int,std::vector<node_ptr>);
    ~CElement();

    void add_node(node_ptr node){m_nodes.push_back(node);}
    void add_id(unsigned int id){m_id=id;}
    int getindex(){return m_id;}
    std::vector<double> get_xcoords();
    std::vector<double> get_ycoords();
    std::vector<unsigned int> get_node_ids();
    boundary_ptr get_boundary(){return m_boundary;}
    node_ptr getnode(unsigned int i) const{return m_nodes.at(i);}
    const std::vector<node_ptr>& get_nodes()const
            {return m_nodes;}
    void show_node_ids(){std::cout<<"\n"; for(auto& iter:m_nodes) std::cout<<"\t"<< iter->get_id();}
    bool has_constrained_face(){if(m_boundary->m_constrained_faces.size()!=0) return true;else false;}
    double get_objval(){return m_objval;}



   // virtual void get_fval   (const double* x,CMaterial*)=0;
   // virtual void set_grad_f (const double* x,const double* g, CMaterial*) =0;
   // virtual int                    get_dof_count()                       const =0;
   // virtual void                   set_node_dof()                       =0;

    virtual void add_boundary(boundary_ptr temp) {m_boundary=temp;}
    virtual void set_boundary_faces() =0;
    virtual     void set_boundary_info()=0;
    virtual void set_constrained_faces()=0;
    virtual void set_boundary_constraints(const double*, double*)=0;
    virtual void set_constraints_boundvals()=0;
    virtual void set_NLP_constraint_matrix(BoundaryFace&)=0;

    virtual std::vector<node_ptr>  get_face(unsigned int)                const =0;    // return node_ptrs of nth face of element.
    virtual Element::element_type           get_element_type()                    const =0;
    virtual int                    faces_count()                         const =0;
    virtual bool                   is_boundary_element()                         const =0;
    virtual int                    elements_shared_by_face(unsigned int) const =0;


    //virtual double                 get_objval(const std::vector<double>&, material_ptr,
   //                                        FEval&)    const =0;
    virtual void                   calculate_fvals(std::vector<double>&) =0;


    boundary_ptr m_boundary;

protected:
    Element::element_type m_type;
    int m_id;
    std::vector<node_ptr> m_nodes;
    double         m_objval;
    shape_ptr m_shape;
    bool is_boundary;



};

#endif // CELEMENT_H
