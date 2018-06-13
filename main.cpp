#include <iostream>
//#include "cmodel.h"
//#include "cfileio.h"
#include<ctime>
#include"cnode.h"
#include"element/celement.h"
#include"element/cps4.h"
#include "cmodel.h"
#include "cfileio.h"
#include "eigen3/Eigen/Core"
#include "NLP/cps4_nlp.h"
#include"cquadratrue.h"
#include <cstdlib>
#include<cmath>
using namespace std;

//from Ipopt library
#include"coin/IpIpoptApplication.hpp"
#include"coin/IpSolveStatistics.hpp"
#include"coin/IpTNLP.hpp"

int main()
{
    clock_t t1;
    t1 =clock();
    t1 =clock() - t1;

    //creating instance of model before reading the input file.
    CModel model;


    //reading the input file.
    {
        CFileio read("input/input.in",&model);
    }

    //create an instance of N.L.P
    SmartPtr<TNLP> problem = new CPS4_NLP(&model);

    // Create an instance of the IpoptApplication
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    //set options.
    app->Options()->SetNumericValue("tol", 1e-9);
    app->Options()->SetStringValue("hessian_approximation","limited-memory");
//    app->Options()->SetStringValue("derivative_test","first-order");
//    app->Options()->SetStringValue("linear_solver","mumps");
//    app->Options()->SetStringValue("dependency_detector","mumps");

    // Intialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      printf("\n\n*** Error during initialization!\n");
      return (int) status;
    }

    // Ask Ipopt to solve the problem
      status = app->OptimizeTNLP(problem);

      if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
      }
      else {
        printf("\n\n*** The problem FAILED!\n");
      }

      // As the SmartPtrs go out of scope, the reference count
      // will be decremented and the objects will automatically
      // be deleted.

      return (int) status;
}
/*
int main()
{


         CModel model;
     CFileio read("input/sample.in",&model);
     auto& temp=  model.get_nodes();
     for(auto&& iter:temp)
     {
         cout<<iter->get_id()<<endl;
     }

     auto& els = model.get_elements();
     for(auto&& it:els)
     {
         cout<<it->getindex()<<endl;
         it->show_node_ids();cout<<"\n";
         cout<<it->is_boundary_element()<<"\n";
     }

     auto & nodes=model.get_nodes();

     for(auto& iter:nodes)
     {
         cout<<iter->get_id()<<"\t";
         if(iter->is_boundary_node()) cout<<"true\n"; else cout<<"false\n";

     }
     model.show_boundary_elements();
     model.boundary_size();
     cout<<"\n\nBOUNDARY\n";
     for(auto& iter:model.m_boundary_elems)
     {
         cout<<iter->getindex()<<"\t"<<iter->has_constrained_face()<<"\n";
     }
//     auto boundary=model.get_boundary_elements();
//     for(auto& iter:boundary)
//     {
//         cout<<iter->getindex();
//     }

}
*/
