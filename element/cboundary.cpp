#include "cboundary.h"
#include <cmath>
BoundaryFace::BoundaryFace()
{
    is_constraint_set=false;
    m_index = -1;

}

BoundaryFace::BoundaryFace(node_ptr n1,node_ptr n2)
{
    node1=n1;
    node2=n2;
    is_constraint_set=false;
    m_index = -1;

    double x2= this->node2->x();
    double x1= node1->x();
    double y2=node2->y();
    double y1= node1->y();

    nx= (x2-x1)/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    ny=(y2-y1)/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));


}


double BoundaryFace::length()
{
    double x2= this->node2->x();
    double x1= node1->x();
    double y2=node2->y();
    double y1= node1->y();

   double length= sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
   return length;
}
