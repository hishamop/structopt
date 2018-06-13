#include "cnode.h"



void CNode::add_incident_elem_id(unsigned int  id)
{
    m_incident_elements.push_back(id);
}

std::vector<unsigned int> CNode::get_shared_elements_id()const
{
    return m_incident_elements;
}

bool CNode::is_boundary_node()
{
    return is_boundary;
}

void CNode::update_dof(const double *x)
{
    if(m_fixity==0)
    {
       m_sdof = &x[(m_id-1)*8];
       m_ddof = &x[(m_id-1)*8+6];
       m_traction =&m_zero[0];
    }
    else
    {
        m_sdof = &x[(m_id-1)*8];
        m_ddof =&m_zero[0];
        m_traction=  &x[(m_id-1)*8+6];
    }
}

std::vector<double> CNode::get_dof()
{
    std::vector<double> dof;
    for(int i=0;i<6;i++)
        dof.push_back(this->m_sdof[i]);
    dof.push_back(m_ddof[0]);
    dof.push_back(m_ddof[1]);
    return dof;


}

void CNode::update_dof(std::vector<double> &x)
{
    m_sdof = &x[(m_id-1)*8];
    m_ddof = &x[(m_id-1)*8+6];
}
