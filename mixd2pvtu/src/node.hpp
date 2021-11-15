#ifndef NODE_H_
#define NODE_H_

#include <stdlib.h>

    class Node
    {
        public:
        Node();
        // Node(Node &&) = default;
        // Node(const Node &) = default;
        // Node &operator=(Node &&) = default;
        // Node &operator=(const Node &) = default;
        ~Node();

        inline void setX(double value) {x = value;};
        inline void setY(double value) {y = value;};
        inline void setZ(double value) {z = value;};

        inline double getX() {return x;};
        inline double getY() {return y;};
        inline double getZ() {return z;};

        inline double getCoord(int i)
        {
            switch(i)
            {
                case 0: return x;
                case 1: return y;
                case 2: return z;
                default: exit(-1);
            }
        };

        private:
        double x;
        double y;
        double z;


    };


#endif // !NODE_H_
