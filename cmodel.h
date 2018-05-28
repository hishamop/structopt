#ifndef CMODEL_H
#define CMODEL_H
#include "cnode.h"
#include "element/celement.h"
#include "element/cps4.h"
#include "material/cmaterial.h"
#include<vector>
#include<list>
#include<memory>
#include<map>


using node_ptr =std::shared_ptr<CNode>;
using elem_ptr = std::shared_ptr<CElement>;
using material_ptr = std::shared_ptr<CMaterial>;
using boundary_ptr = std::shared_ptr<CBoundary>;
class CModel
{

    friend class CFileio;
    friend class CNlp;

public:
    CModel();

    //INTERFACE FUNCTIONS
    void add_material(std::string,material_ptr);
    void set_dof();



    void                        add_node_ptr(node_ptr t)
    {
        m_nodes.push_back(t);
    }
    void                        add_node(node_ptr t)
    {
        m_nodes.push_back(t);
    }
    void                        add_element(elem_ptr t)
    {
        m_elements.push_back(t);
    }

    const std::vector<elem_ptr>& get_element_vector() const
    {
        return m_elements;
    }

    node_ptr                    get_node_ptr(int id) const{
        return m_nodes.at(id);
    }

    void add_boundary_element(elem_ptr t)
    {
        m_boundary_elems.push_back(t);
    }

    elem_ptr  get_elem_ptr(unsigned int i) const{return m_elements.at(i);}

    std::vector<node_ptr>& get_nodes(){return m_nodes;}
    std::vector<elem_ptr>& get_elements(){return m_elements;}
    std::vector<elem_ptr>& get_boundary_elements(){return m_boundary_elems;}
    void show_boundary_elements(){if( m_boundary_elems.size()==0) std::cout<<"\n Boundary is not set"; std::cout<<"\n";for(auto& iter:m_boundary_elems) std::cout<<"\t"<<iter->getindex();std::cout<<"\n";}
    void show_nset(){}
    void show_elsets(){}
    void boundary_size(){std::cout<<"\n"<<m_boundary_elems.size();}

public:
    std::vector<node_ptr>                            m_nodes;
    std::vector<elem_ptr>                            m_elements;
    std::vector<elem_ptr>                            m_boundary_elems;
//    std::map<std::string,std::vector<unsigned int>>  m_nodeset_map;
    std::map<std::string,std::vector<node_ptr>>  m_nodeset_map;

    //std::map<std::string,std::vector<unsigned int>>  m_elset_map;
    std::map<std::string,std::vector<elem_ptr>>  m_elset_map;
    material_ptr                                     m_material;
   // std::map<unsigned int, boundary_ptr>             m_boundary_map;
  //  std::vector<boundary_ptr>                        m_boundary;


    CModel (CModel&);    //Compiler defined copy const. not allowed. So defined as private.

};

#endif // CMODEL_H



