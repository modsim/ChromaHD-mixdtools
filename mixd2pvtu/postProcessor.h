#ifndef POSTPROCESSOR_H_
#define POSTPROCESSOR_H_

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
/* #include "LDEMVTKXMLPUnstructuredDataWriter.h" */
#include <vtkSmartPointer.h>
#include "settings.h"
#include "tet.h"

class postProcessor
{
    private:
        //VARIABLES
        inputSettings*  settings;   // a local pointer for the settings
        tetMesh*        mesh;       // a local pointer for the mesh
        double          partialminT;       // min value of the Temperature field
        double          partialmaxT;

        //METHODS
        void vtkVisualization(int irec);
    protected:

    public:
        void postProcessorControl(inputSettings*, tetMesh*, int);
};


#endif /* POSTPROCESSOR_H_ */
