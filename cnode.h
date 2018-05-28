#ifndef CNODE_H
#define CNODE_H
#include"eigen3/Eigen/Core"
#include<vector>



class CNode
{
public:
    CNode():m_zero{0.0,0.0,0.0}{m_id=0;m_x=0.0;m_y=0.0;is_boundary=false; for(int i=0;i<3;i++);}
    CNode(int id,double x, double y):m_id(id),m_x(x),m_y(y),m_zero{0.0,0.0,0.0}{is_boundary=false;}
    unsigned int get_id()    {return m_id;}
    void update_stress_dof(const double*);
    void update_disp_dof(const double*);
    void update_traction(const double*);
    void update_dof(const double*);
    void update_grad_dof();
    std::vector<double>& get_gradients(){return m_grad_dof;}
    double x(){return m_x;}
    double y(){return m_y;}

    std::vector<unsigned int> get_shared_elements_id()const ;
    bool is_boundary_node();
    bool is_constrined(){if(m_fixity==0) return false;else return true;}
    void on_boundary(){is_boundary=true;}
    void add_incident_elem_id(unsigned int);
    void set_fixity(int f){m_fixity=f;}
    const double* get_stress_dof()const{return m_sdof;}
    const double* get_displacement_dof()const{return m_ddof;}

    double sdof(unsigned int i){return m_sdof[i];}
    double ddof(unsigned int i){return m_ddof[i];}

    //stress at node
    double stress_xx();
    double stress_yy();
    double stress_xy();


private:



    unsigned int m_id;
    double m_x;
    double m_y;
    std::vector<unsigned int> m_incident_elements;
    unsigned int m_fixity;
    unsigned int   m_Ndof;
    std::vector<double> m_grad_dof;
    const double* m_traction;
    bool is_boundary;
    const double* m_sdof;
    const double* m_ddof;
    std::vector<double> m_stsdof;
    std::vector<double>m_dispdof;
    const double m_zero[3];

};

#endif // CNODE_H
