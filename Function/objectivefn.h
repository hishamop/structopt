#ifndef OBJECTIVEFN_H
#define OBJECTIVEFN_H
#include "cmodel.h"

class ObjectiveFn
{
public:
    ObjectiveFn(CModel*);
    void get_objval(std::vector<double>);

private:
    double m_objval;
    CModel* m_model;
    int m_nvar;
};

#endif // OBJECTIVEFN_H
