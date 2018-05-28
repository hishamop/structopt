#include "cfileio.h"
#include "utility/auxil.h"
#include <algorithm>
#include <iomanip>
#include <assert.h>
#include"element/cboundary.h"


bool custom_compare(elem_ptr obj, unsigned int id)  { return obj->getindex() < id; }


CFileio::CFileio(std::string szFile, CModel *model)
{
    {
        m_model = model;
        m_File.open(szFile,std::ios::in);
        if(!m_File)
        std::cout<<"Coudn't open the file";

        //read the file line by line
        std::string s;
        while(std::getline(m_File,s))
        {
            //convert all letters into upper.
            std::string upper(s);
            std::transform(upper.begin(),upper.end(),upper.begin(),toupper);
            upper.erase(std::remove_if(upper.begin(), upper.end(), isspace), upper.end());
            if(upper.find("*NODE") != std::string::npos)
            {
                read_nodes();
            }

            if(upper.find("*ELEMENT") != std::string::npos)
            {
                read_elems(upper);
            }

            if(upper.find("*NSET") != std::string::npos)
            {
                read_nodeset(upper);

            }

            if(upper.find("*ELSET") != std::string::npos)
            {
                read_elset(upper);

            }
           if(upper.find("*STEP") != std::string::npos)
            {
                set_boundary();
            }

            if(upper.find("*PLOAD") != std::string::npos)
            {
                read_pload();
            }


            if(upper.find("*BOUNDARY") != std::string::npos)
            {
                read_boundary();
            }
        }

        for(auto& iter:m_model->m_boundary_elems)
        {
            iter->set_constrained_faces();
        }

    }
}



//prive class functions
void CFileio::read_nodes()
{
    char c;                         //temp.variable for recieving ','.
    int id;
    float x,y,z;
    //The program will fail if the input file contains blank lines.
    while(m_File.peek() != '*' and m_File.peek() != EOF )
    {
        if(m_File.peek() == '\n')              //if the line is blank skip to next line.
        {
            std::string skipline;
            std::getline(m_File,skipline); continue;
        }

        if(m_File>>id>>c>>x>>c>>y)              //input file consist only x,y coordinates
        {
            //node temp(id,x,y);
            node_ptr temp_ptr = std::make_shared<CNode>(id,x,y);
            m_model->add_node(temp_ptr);
            std::string dummy;
            std::getline(m_File,dummy);
            continue;
        }

//        if(m_File>>id>>c>>x>>c>>y>>c>>z)           //input file consist x,y,z coordinates.
//        {
//            //node temp(id,x,y,z);
//            node_ptr temp_ptr =std::make_shared<CNode>(id,x,y,z);
//            m_model->add_node(temp_ptr);
//         //   m_model->add_node(temp);
//            std::string dummy;
//            std::getline(m_File,dummy);
//            continue;
//        }
    }
}




void CFileio::read_nodeset(std::string line)
{
    std::string nodeset_name = parse_label(line,"NSET");
//    nset_ptr temp = std::make_shared <nodeset> (nodeset_name);
    std::vector<unsigned int> id_s;
    while(m_File.peek() != '*' and m_File.peek() != EOF )
    {
        if(m_File.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_File,dummy);continue;        // skip the blank line.
        }
        std::string ss; int id;
        std::getline(m_File,ss);
        std::stringstream str(ss);


        std::vector<node_ptr> node_list;
        while(str>>id)
        {
            id_s.push_back(id);
            node_list.push_back(m_model->m_nodes.at(id-1));
//            temp->add_nset_node(m_model->get_node_ptr(id));
            if(str.peek() == ',' or str.peek() == ' ')
            {
                str.ignore();
            }

        }

        m_model->m_nodeset_map[nodeset_name] = node_list;

    }
 //   m_model->add_nodeset(temp);
}


void CFileio::read_elset(std::string line)
{
    std::string elset_name = parse_label(line,"ELSET");
//    elset_ptr temp = std::make_shared<elset>(elset_name);
    while(m_File.peek() != '*' and m_File.peek() != EOF )
    {
        if(m_File.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_File,dummy);continue;        // skip the blank line.
        }
        std::string ss; int id;
        std::getline(m_File,ss);
        std::transform(ss.begin(),ss.end(),ss.begin(),toupper);
        std::stringstream str(ss);
        std::vector<unsigned int> element_ids;
        std::vector<elem_ptr> element_ptr_list;
        while(str>>id)
        {
            element_ids.push_back(id);

            element_ptr_list.push_back(m_model->m_elements.at(id-1));
            if(str.peek() == ',' or str.peek() == ' ')
            {
                str.ignore();
            }

        }
        m_model->m_elset_map[elset_name] = element_ptr_list;
    }

    //m_model->add_elset(temp);
}

//read elements from the input fle and store it in m_element.
void CFileio::read_elems(std::string upper)
{
    int id; char c;
    //below loop for quad4 element
    if(upper.find("CPS4") !=std::string::npos)
    {
        while(m_File.peek() != '*' and m_File.peek() != EOF)
        {
            if(m_File.peek() == '\n')
            {
                std::string dummy;
                std::getline(m_File,dummy); continue;    // skip the blank line.
            }
            int n[4];
            if(m_File>> id >> c >> n[0]>> c >> n[1]>> c >> n[2]>> c >>n[3])
            {
                elem_ptr temp_elem =  std::make_shared<CPS4>(id);
                for(int i=0; i<4; i++)
                {
                    temp_elem->add_node(m_model->m_nodes.at(n[i]-1));
                    m_model->m_nodes.at(n[i]-1)->add_incident_elem_id(id);
                }
                m_model->add_element(temp_elem);
                std::string dummy;
                std::getline(m_File,dummy);
                continue;  // dummy string to collect junk data is not used
            }
        }
    }
}


void CFileio::read_boundary()
{
    while(m_File.peek() != '*' and m_File.peek() != EOF )
    {
        if(m_File.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_File,dummy);continue;        // skip the blank line.
        }
        std::string ss;
        std::getline(m_File,ss);


        // Remove whitespace.
        ss.erase(std::remove_if(ss.begin(), ss.end(), isspace), ss.end());
        std::transform(ss.begin(),ss.end(),ss.begin(),toupper);

        std::string nset_name; unsigned int f1,f2;
        {
            char comma =','; char c;
            std::stringstream ss_stream(ss);
            {
                int flag =0;
                while(ss_stream>>c)
                {
                    if(c==comma)      break;
                }
                ss_stream>>f1>>c>>f2;
            }

            size_t first_comma_index;
             first_comma_index = ss.find(",");
            {
                std::string::iterator
                        beg =ss.begin(),
                        end =ss.begin()+first_comma_index;
                nset_name= std::string(beg,end);
            }

        }

        for(auto& iter:m_model->m_nodeset_map[nset_name])
        {
            iter->set_fixity(3);

        }

    }


}


/*void CFileio::read_nloads()
{
    while(m_File.peek() != '*' and m_File.peek() != EOF )
    {
        if(m_File.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_File,dummy);continue;        // skip the blank line.
        }
        std::string ss;
        std::getline(m_File,ss);
        std::stringstream str_id(ss);
        std::stringstream str_name(ss);

        char c;
        std::string elset_name;
        int elem_id;
        unsigned int side;
        unsigned int orientation;
        double tr_n1;   // Traction at first node
        double tr_n2;   //traction at last node.

        if(str_name>>elset_name)
        {
            str_name>>c>>side>>c>>orientation>>c>>tr_n1>>c>>tr_n2;

        }

        if(str_id>>elem_id)
        {
            str_id>>c>>side>>c>>orientation>>c>>tr_n1>>c>>tr_n2;
                std::cerr<<"Boundary edge not found";

        }

    }
}*/





 void CFileio::set_boundary()
{
    static bool flag = true;
    if(flag)
    {
        for(auto& iter: m_model->m_elements)
        {
            if(iter->is_boundary_element())
            {
              //  m_model->m_boundary_elems.push_back(iter);
                    boundary_ptr temp =std::make_shared<CPS4_boundary>();

                  //  m_model->m_boundary.push_back(temp);
                   // temp->set_boundary_faces();

                    iter->add_boundary(temp);
                    iter->set_boundary_info();
                    m_model->m_boundary_elems.push_back(iter);
                   // iter->set_boundary_faces();
                    //m_model->add_boundary_element(temp);
                    //m_model->m_boundary_map[iter->getindex()] = temp;
                    //break;)
            }
        }

        flag = false;
    }
}



/*void CFileio::read_pload()
{

}
*/



void CFileio::read_pload()
{
    while(m_File.peek() != '*' and m_File.peek() != EOF )
    {
        if(m_File.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_File,dummy);continue;        // skip the blank line.
        }

        std::string ss;   // line containg elset,face,pressure
        std::getline(m_File,ss);

        // Remove whitespace.
        ss.erase(std::remove_if(ss.begin(), ss.end(), isspace), ss.end());
        std::transform(ss.begin(),ss.end(),ss.begin(),toupper);

        std::string elset_name, szface; double pressure;
        {
            char comma =','; char c;
            std::stringstream ss_stream(ss);
            {
                int flag =0;
                while(ss_stream>>c)
                {

                    if(c==comma)
                       {
                           flag+=1;
                           if(flag ==2)
                               break;
                    }

                }
                ss_stream>>pressure;        //get pressure load.
            }

            size_t first_comma_index, second_comma_index;
             first_comma_index = ss.find(",");
            second_comma_index = ss.find(",",first_comma_index);
            {
                std::string::iterator
                        beg =ss.begin(),
                        end =ss.begin()+first_comma_index;
                elset_name= std::string(beg,end);
            }
            {
               std::string::iterator
                        beg =ss.begin() +first_comma_index+1,
                        end =ss.begin() +second_comma_index+first_comma_index-1;
                szface =std::string(beg,end);
            }   //get element set and load orientation

        }

        unsigned int nface=0;
        if(szface=="P2")
            nface =1;
        else if(szface=="P3")
            nface =2;
        if(szface =="P4")
            nface =3;


        for(auto& iter:m_model->m_elset_map[elset_name])
        {
            iter->get_boundary()->m_pload[nface]=pressure;

        }


        //Assign load and orientation to concerned elements.
      /*  {
            assert(m_model->m_elset_map.count(elset_name)>0);


            for(auto  &iter: m_model->m_elset_map[elset_name])
            {

                load_ptr temp = std::make_shared<Dload>(nface,pressure);
                auto boundary = m_model->m_boundary_map[iter];

                boundary->set_load(temp);
                boundary->set_normal_pressure_on_face(nface,pressure);
//               m_model->m_boundary_map[iter]->set_normal_pressure_on_face(nface,pressure);
            }

        }*/
    }
}















