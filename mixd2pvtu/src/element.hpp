#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <string>
#include <vtk/vtkCellType.h>
#include <iostream>

class Element{
    public:
        virtual void setConn  (int i, int value) {conn[i] = value;};
        virtual void setLConn (int i, int value) {lConn[i] = value;};

        virtual int getConn  (int index) {return conn[index];};
        virtual int getLConn (int index) {return lConn[index];};

        int * conn;
        int * lConn;
};

class Tetrahedron
    :public Element
    {

    public:

    Tetrahedron(){
        conn = new int[4];
        lConn = new int[4];
        };
    ~Tetrahedron(){
        delete[] conn;
        delete[] lConn;
    };

};

class Triangle
    :public Element
    {

    public:

    Triangle(){
        conn = new int[3];
        lConn = new int[3];
        };

    ~Triangle(){
        delete[] conn;
        delete[] lConn;
    };

};

#endif
