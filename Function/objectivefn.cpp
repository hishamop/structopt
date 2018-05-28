#include "objectivefn.h"
#include"feval.h"
#include "material/cmaterial.h"
ObjectiveFn::ObjectiveFn(CModel* model)
{
    m_objval=0;
    m_model=model;
    m_nvar=m_model->m_nodes.size()*8;
}




void ObjectiveFn::get_objval()
{

    FEval preval(CPS4);
    STEEL steel;
    for(auto& iter:m_model->m_elements)
    {

        m_objval+=iter->get_fval();
    }

}
