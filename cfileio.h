#ifndef CFILEIO_H
#define CFILEIO_H
#include <iostream>
#include <fstream>
#include "cmodel.h"
#include "cnode.h"
class CModel;
class CFileio
{
public:
    CFileio(std::string , CModel*);

private:
    std::ifstream m_File;
    CModel* m_model;
    void read_nodes();
    void read_elems(std::string);
    void read_nodeset(std::string);
    void read_elset(std::string);
   // void read_nloads();
   // void read_dload();
    void read_pload();
   // void read_disp();
   // void read_fixity();
    void set_boundary();
    void read_boundary();
};

#endif // CFILEIO_H
