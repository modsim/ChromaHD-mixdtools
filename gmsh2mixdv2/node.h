#ifndef NODES_H
#define NODES_H

#include<string>
#include<sstream>

class Node
{
    public:
        double x;
        double y;
        double z;

        Node(double x_, double y_, double z_)
        {
            x = x_;
            y = y_;
            z = z_;
        }

        Node(std::string line)
        {
            std::istringstream iss(line);
            size_t id;
            iss >> id >> x >> y >> z;
        }

};

#endif
