from abaqus import *
from abaqusConstants import *

def partition_matching_holes():
    # 1. Connect to your model and part 
    # (Update 'Model-1' and 'Part-1' to match your actual Abaqus database names)
    my_model = mdb.models['Model-1']
    my_part = my_model.parts['Part-1']
    
    target_radius = 5.0
    tolerance = 0.01 
    
    # 2. Iterate through all faces to find cylindrical holes
    for face in my_part.faces:
        try:
            # getSurfaceCylinder() returns (PointOnAxis, AxisDirection, Radius).
            # It intentionally raises an exception if the face is NOT a cylinder.
            cylinder_data = face.getSurfaceCylinder()
            radius = cylinder_data[2]
            
            # 3. Check if the radius matches our target
            if abs(radius - target_radius) < tolerance:
                center_point = cylinder_data[0]
                print("Matching hole found at: {}".format(center_point))
                
                # 4. Create Datum Planes at the hole's center
                # This creates standard principal planes passing through the hole's center
                plane_xy = my_part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=center_point[2])
                plane_xz = my_part.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=center_point[1])
                
                # 5. Partition all cells in the part using the new datum planes
                # We use a try/except here because partitioning a cell that is already 
                # fully split by a previous loop iteration can sometimes throw a warning.
                try:
                    my_part.PartitionCellByDatumPlane(cells=my_part.cells, datumPlane=my_part.datums[plane_xy.id])
                except:
                    pass
                    
                try:
                    my_part.PartitionCellByDatumPlane(cells=my_part.cells, datumPlane=my_part.datums[plane_xz.id])
                except:
                    pass
                    
        except:
            # The face is not a cylinder; skip it and move to the next face
            pass
            
    print("Partitioning complete.")

# Execute the function
partition_matching_holes()