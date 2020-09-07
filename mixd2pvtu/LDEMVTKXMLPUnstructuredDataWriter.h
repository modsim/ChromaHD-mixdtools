#include "vtkXMLPUnstructuredGridWriter.h"

#include "vtkObjectFactory.h"
#include "vtkSetGet.h"

class LDEMVTKXMLPUnstructuredDataWriter : public vtkXMLPUnstructuredGridWriter
{
public:

    inline static LDEMVTKXMLPUnstructuredDataWriter* New();

    void WritePPieceAttributes(int index)
    {
        std::string DATAFNAME = "test-case_";
        DATAFNAME += std::to_string(index);
        /* DATAFNAME += "/spheredata_"; */
        /* DATAFNAME += std::to_string(index); */
        DATAFNAME += ".vtu";

        this->WriteStringAttribute("Source", DATAFNAME.c_str());
    }

};
vtkStandardNewMacro(LDEMVTKXMLPUnstructuredDataWriter);
