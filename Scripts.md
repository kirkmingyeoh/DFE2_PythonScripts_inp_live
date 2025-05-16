# Scripts

This document lists the Python scripts available in the repository and describes the corresponding Direct FE2 model. The scripts are generally named in the format 'DFE2_Dimension_MacroEle-FurtherInfo_RVEEle-FurtherInfo_FurtherModelInfo.py'. 

Unless otherwise stated, these scripts assume the following:  
Macroscale - all macroscale elements are homogenised with the same RVE and same number of integration points  
RVE - rectangular (2D) or cuboidal (3D) geometry; perfectly periodic mesh  

| Name | Dimension | Macro element | RVE element | Additional details |
| :-----: | :-----: | :-----: | :-----: | :-----: |
| DFE2_2D_Quad_QuadTri.py | 2D | Quadrilateral | Quadrilateral or triangle | - |
| DFE2_3D_Hex_HexTetWedge.py | 3D | Hexahedral | Hexahedral, wedge or tetrahedral | - |

