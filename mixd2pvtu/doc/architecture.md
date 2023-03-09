---
pagestyle: empty
geometry:
    - top=0mm
    - left=0mm
    - bottom=0mm 
    - right=0mm
---

\centering

# Current/Old

```{.mermaid format=pdf}
flowchart TB
    set[Settings.Prepare]
    meshprep[Mesh.Prepare]
    meshproc[Mesh.ProcessData]
    post[PostProcess]
    
    set --> meshprep --> meshproc --> post
    
    subgraph SETTINGS
    direction TB
    A[Read File] --> B[Read Arguments] --> C[Print Settings]
    end
    
    set -.-> SETTINGS
    
    subgraph MESH
    direction TB
    rmf[Read mesh] --> gts[Get Timesteps] --> loc[Localize coords]
    end
    
    meshprep -.-> MESH
    
    subgraph DATA
    direction TB
    read[Read Data] --> locd[Localize Data]
    end
    
    meshproc -.-> DATA
```

# Classes

```{.mermaid format=pdf}
classDiagram
    class Config {
        +Title = "vis"
        +Outpath = "output"
        +minf = "minf"
        +mien = "mien"
        +mxyz = "mxyz"
        +datafiles[] =  
        +ndf = 0
        +nrec = 0
        +nrecoffset = 0
        +nrecstride = 1
        +meshtype = ST or SD
        +dt = 1
        +dtFile = ""
        +ReadFile()
        +ReadCmdline()
    }
    Element <|-- Tet
    Element <|-- Tri
    Mesh --> Element : contains
    Mesh --> Node : contains
    class Element {
    }
    class Mesh {
        +Data
        +Nodes
        +Elements
        +read()
        +localizeNodes()
        +write()
    }
    class Tet {
        +VTKElemType = VTK_TETRA
        +nen = 4
        +nef = 4
    }
    class Tri {
        +VTKElemType = VTK_TRIANGLE
        +nen = 3
        +nef = 2 or 3?
    }

```

