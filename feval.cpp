#include "feval.h"
#include"math.h"



void FEval::compute()
{
    switch (m_type) {
    case Element::CPS4:
        cps4_bsmat();
        break;
    default:
        break;
    }
}

void FEval::cps4_bsmat()
{
    for(auto sval:m_quad.get_nodes())
    {
        for (auto tval:m_quad.get_nodes())
        {
            Eigen::MatrixXd Bs =Eigen::MatrixXd::Zero(3,24);
            Bs(0,0)=(3*tval*(sval - 1)*(pow(sval,2) + sval + 5*pow(tval,2) - 5))/8;
            Bs(0,1)= 3*tval*pow(sval -1 , 2)*(sval + 1)/8;
            Bs(0,2) = (sval - 1)*(3*pow(sval , 2)*tval - pow(sval , 2) + 3*sval*tval - sval + 15*pow(tval,3) - 3*pow(tval,2) - 15*tval + 3)/8;
            Bs(0,4) = (sval - 1)*(tval - 1)*(5*pow(tval,2) + 2*tval - 1)/8;
            Bs(0,5) = (3*tval - 1)*pow(sval -1 , 2)*(sval + 1)/8;
            Bs(0,6) = 3*tval*(sval + 1)*(- pow(sval , 2) + sval - 5*pow(tval,2) + 5)/8;
            Bs(0,7) = 3*tval*(sval - 1)*pow(sval +1 , 2)/8;
            Bs(0,8) = -(sval + 1)*(3*pow(sval , 2)*tval - pow(sval , 2) - 3*sval*tval + sval + 15*pow(tval,3) - 3*pow(tval,2) - 15*tval + 3)/8;
            Bs(0,10) = -(sval + 1)*(tval - 1)*(5*pow(tval,2) + 2*tval - 1)/8;
            Bs(0,11) = (3*tval - 1)*(sval - 1)*pow(sval +1 , 2)/8;
            Bs(0,12) = 3*tval*(sval + 1)*(pow(sval , 2) - sval + 5*pow(tval,2) - 5)/8;
            Bs(0,13) = -3*tval*(sval - 1)*pow(sval +1 , 2)/8;
            Bs(0,14) = (sval + 1)*(-3*pow(sval , 2)*tval - pow(sval , 2) + 3*sval*tval + sval - 15*pow(tval,3) - 3*pow(tval,2) + 15*tval + 3)/8;
            Bs(0,16) = (sval + 1)*(tval + 1)*(5*pow(tval,2) - 2*tval - 1)/8;
            Bs(0,17) = (3*tval + 1)*(sval - 1)*pow(sval +1 , 2)/8;
            Bs(0,18) = -(3*tval*(sval - 1)*(pow(sval , 2) + sval + 5*pow(tval,2) - 5))/8;
            Bs(0,19) = -3*tval*pow(sval -1 , 2)*(sval + 1)/8;
            Bs(0,20) = (sval - 1)*(3*pow(sval , 2)*tval + pow(sval , 2) + 3*sval*tval + sval + 15*pow(tval,3) + 3*pow(tval,2) - 15*tval - 3)/8;
            Bs(0,22) = (sval - 1)*(tval + 1)*(- 5*pow(tval,2) + 2*tval + 1)/8;
            Bs(0,23) = (3*tval + 1)*pow(sval -1 , 2)*(sval + 1)/8;

            Bs(1,0) = (3*sval*(tval - 1)*(5*pow(sval , 2) + pow(tval,2) + tval - 5))/8;
            Bs(1,1) = (tval - 1)*(15*pow(sval,3) - 3*pow(sval , 2) + 3*sval*pow(tval,2) + 3*sval*tval - 15*sval - pow(tval,2) - tval + 3)/8;
            Bs(1,2) = 3*sval*pow((tval - 1) , 2)*(tval + 1)/8;
            Bs(1,3) = (sval - 1)*(tval - 1)*(5*pow(sval , 2) + 2*sval - 1)/8;
            Bs(1,5) = (3*sval - 1)*pow((tval - 1) , 2)*(tval + 1)/8;
            Bs(1,6) = -3*sval*(tval - 1)*(5*pow(sval , 2) + pow(tval,2) + tval - 5)/8;
            Bs(1,7) = (tval - 1)*(15*pow(sval,3) + 3*pow(sval , 2) + 3*sval*pow(tval,2) + 3*sval*tval - 15*sval + pow(tval,2) + tval - 3)/8;
            Bs(1,8)= -3*sval*pow((tval - 1) , 2)*(tval + 1)/8;
            Bs(1,9)= ((sval + 1)*(tval - 1)*(- 5*pow(sval , 2) + 2*sval + 1))/8;
            Bs(1,11) = (3*sval + 1)*pow((tval - 1) , 2)*(tval + 1)/8;
            Bs(1,12) = 3*sval*(tval + 1)*(5*pow(sval , 2) + pow(tval,2) - tval - 5)/8;
            Bs(1,13) = (tval + 1)*(-15*pow(sval,3) - 3*pow(sval , 2) - 3*sval*pow(tval,2) + 3*sval*tval + 15*sval - pow(tval,2) + tval + 3)/8;
            Bs(1,14)= -3*sval*(tval - 1)*pow((tval + 1) , 2)/8;
            Bs(1,15) = (sval + 1)*(tval + 1)*(5*pow(sval , 2) - 2*sval - 1)/8;
            Bs(1,17) = (3*sval + 1)*(tval - 1)*pow((tval + 1) , 2)/8;
            Bs(1,18) = (3*sval*(tval + 1)*(-5*pow(sval , 2) - pow(tval,2) + tval + 5))/8;
            Bs(1,19) = (tval + 1)*(-15*pow(sval,3) + 3*pow(sval , 2) - 3*sval*pow(tval,2) + 3*sval*tval + 15*sval + pow(tval,2) - tval - 3)/8;
            Bs(1,20) = 3*sval*(tval - 1)*pow((tval + 1) , 2)/8;
            Bs(1,21) = -(sval - 1)*(tval + 1)*(5*pow(sval , 2) + 2*sval - 1)/8;
            Bs(1,23) = (3*sval - 1)*(tval - 1)*pow((tval + 1) , 2)/8;

            Bs(2,0) = -(15*pow(sval , 4) + 18*pow(sval , 2)*pow(tval,2) - 36*pow(sval , 2) + 15*pow(tval , 4) - 36*pow(tval,2) + 24)/32;
            Bs(2,1) = -(3*sval + 1)*(sval - 1)*(5*pow(sval , 2) + 2*sval + 6*pow(tval,2) - 9)/32;
            Bs(2,2) = -(3*tval + 1)*(tval - 1)*(6*pow(sval , 2) + 5*pow(tval,2) + 2*tval - 9)/32;
            Bs(2,3) = -(5*sval + 1)*pow(sval -1 , 2)*(sval + 1)/32;
            Bs(2,4) = -(5*tval + 1)*pow((tval - 1) , 2)*(tval + 1)/32;
            Bs(2,5) = -(3*sval + 1)*(3*tval + 1)*(sval - 1)*(tval - 1)/16;
            Bs(2,6) = (15*pow(sval , 4) + 18*pow(sval , 2)*pow(tval,2) - 36*pow(sval , 2) + 15*pow(tval , 4) - 36*pow(tval,2) + 24)/32;
            Bs(2,7) = (3*sval - 1)*(sval + 1)*(- 5*pow(sval , 2) + 2*sval - 6*pow(tval,2) + 9)/32;
            Bs(2,8) = (3*tval + 1)*(tval - 1)*(6*pow(sval , 2) + 5*pow(tval,2) + 2*tval - 9)/32;
            Bs(2,9) = (5*sval - 1)*(sval - 1)*pow(sval +1 , 2)/32;
            Bs(2,10) = (5*tval + 1)*pow((tval - 1) , 2)*(tval + 1)/32;
            Bs(2,11) = -(3*sval - 1)*(3*tval + 1)*(sval + 1)*(tval - 1)/16;
            Bs(2,12) = -(15*pow(sval , 4) + 18*pow(sval , 2)*pow(tval,2) - 36*pow(sval , 2) + 15*pow(tval , 4) - 36*pow(tval,2) + 24)/32;
            Bs(2,13) = (3*sval - 1)*(sval + 1)*(5*pow(sval , 2) - 2*sval + 6*pow(tval,2) - 9)/32;
            Bs(2,14) = (3*tval - 1)*(tval + 1)*(6*pow(sval , 2) + 5*pow(tval,2) - 2*tval - 9)/32;
            Bs(2,15) = -(5*sval - 1)*(sval - 1)*pow(sval +1 , 2)/32;
            Bs(2,16) = -(5*tval - 1)*(tval - 1)*pow((tval + 1) , 2)/32;
            Bs(2,17) = -(3*sval - 1)*(3*tval - 1)*(sval + 1)*(tval + 1)/16;
            Bs(2,18) = (15*pow(sval , 4) + 18*pow(sval , 2)*pow(tval,2) - 36*pow(sval , 2) + 15*pow(tval , 4) - 36*pow(tval,2) + 24)/32;
            Bs(2,19) = (3*sval + 1)*(sval - 1)*(5*pow(sval , 2) + 2*sval + 6*pow(tval,2) - 9)/32;
            Bs(2,20) = (3*tval - 1)*(tval + 1)*(-6*pow(sval , 2) - 5*pow(tval,2) + 2*tval + 9)/32;
            Bs(2,21) = (5*sval + 1)*pow(sval -1 , 2)*(sval + 1)/32;
            Bs(2,22) = (5*tval - 1)*(tval - 1)*pow((tval + 1) , 2)/32;
            Bs(2,23) = -(3*sval + 1)*(3*tval - 1)*(sval - 1)*(tval + 1)/16;
            m_bsmat.push_back(Bs);
        }
    }
}

void FEval::cps4_bdmat()
{

}
