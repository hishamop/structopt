#include "celement.h"




CElement::~CElement()
{

}

CElement::CElement(int id, std::vector<node_ptr> nlist)
{
    m_id = id;
    m_nodes = nlist;
    is_boundary =false;
    m_shape=nullptr;
}



std::vector<double> CElement::get_xcoords()
{
    std::vector<double> xcoords;
    xcoords.clear();
    xcoords.reserve(4);
    for(auto iter:m_nodes)
    {double xp=iter->x();
        xcoords.push_back(xp);}
    return xcoords;
}

std::vector<double> CElement::get_ycoords()
{
    std::vector<double> ycoords;
    ycoords.clear();
    ycoords.reserve(4);
    for(auto&iter:m_nodes)
        ycoords.push_back(iter->y());
    return ycoords;
}

std::vector<unsigned int>CElement::get_node_ids()
{
    std::vector<unsigned int> Ids;
    Ids.clear();Ids.reserve(4);
    for(auto&iter:m_nodes)
        Ids.push_back(iter->get_id());
    return Ids;
}
