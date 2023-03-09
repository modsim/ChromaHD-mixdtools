#ifndef PHYSICALGROUP_H
#define PHYSICALGROUP_H

#include <string>

class PhysicalGroup
{
    public:
        int dim;
        int tag;
        int matID;
        std::string name;

        PhysicalGroup(int dim_, int tag_, std::string name_, int matID_)
        {
            dim = dim_;
            tag = tag_;
            name = name_;
            matID = matID_;

        }

        ~PhysicalGroup();


};

#endif
