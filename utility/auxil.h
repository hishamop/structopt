#ifndef AUXIL_H
#define AUXIL_H
#include<string>
#include <sstream>
#include<iterator>
#include<algorithm>
std::string parse_label(std::string line,std::string label)
{
    // Remove whitespace.
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

    //Get index of start of "label="
    size_t label_index = line.find(label + "=");

    //return string = str
    std::string str;
    if (label_index != std::string::npos)
      {
        // Location of the first comma following "label="
        size_t comma_index = line.find(",", label_index);
        if(comma_index!=std::string::npos)
        {
            std::string::iterator
                    beg = line.begin() + label.size() + 1 + label_index,
                    end = (comma_index == std::string::npos) ? line.end() : line.begin() + comma_index;
            str = std::string(beg, end);
        }
        else
        {
            std::string::iterator
                    beg = line.begin() + label.size() + 1 + label_index,
                    end = line.end();
            str = std::string(beg, end);
        }

      }
    return str;
}


#endif // AUXIL_H
