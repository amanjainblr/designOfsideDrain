# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# area_calculation.py
# Created on: 2020-12-05 11:37:05.00000,@AMAN_BAGRECHA
#   
# Description: 
# ---------------------------------------------------------------------------

# Set the necessary product code
# import arcinfo


# Import arcpy module
import arcpy
import os

# Local variables:
design_road = r"D:\6th_sem\final_year\arcGIS\arcpy_currentWS.gdb\design_road_Project"
intersecting_parcels = design_road
watershed_raster = r"D:\6th_sem\final_year\arcGIS\arcpy_currentWS.gdb\Watersh_Flow1"
watershed_polygon = "D:\\6th_sem\\final_year\\arcGIS\\currentWS_v2.gdb\\raster_poly_simply_multipart"
currentWS_v2_gdb = "D:\\6th_sem\\final_year\\arcGIS\\currentWS_v2.gdb"

#parcels_divided_by_design_road = "D:\\6th_sem\\final_year\\arcGIS\\currentWS_v2.gdb\\karnataka_highway_Clip_Featu"
#parcels_with_elevation_data = parcels_divided_by_design_road
#Input_Join_Field = "FID_demo"
#Zone_field = "FID_demo"
SRTM_Digital_Elevation_Data_30m_tif = r"D:\6th_sem\final_year\arcGIS\data\dem bangalore\SRTM Digital Elevation Data 30m.tif"
Statistics_type = "MEAN"
#zonal_stats_elevation = "D:\\6th_sem\\final_year\\arcGIS\\currentWS_v2.gdb\\ZonalSt_karnata1"
#Output_Join_Field = "OBJECTID"
#Join_Fields = "MEAN"

# Process: Raster to Polygon

tempEnvironment6 = arcpy.env.scratchWorkspace
arcpy.env.scratchWorkspace = "D:\6th_sem\final_year\arcGIS\area_calculation.gdb"

# overwriting
arcpy.env.overwriteOutput = True

arcpy.env.workspace = r"D:\6th_sem\final_year\arcGIS\area_calculation.gdb"

wkspace= r"D:\6th_sem\final_year\arcGIS\arcpy_output"


# Name: RasterToPolygon_Ex_02.py
# Description: Converts a raster dataset to polygon features.
# Requirements: None

# Import system modules
import arcpy
from arcpy import env



# Set local variables
inRaster = watershed_raster
outPolygons = os.path.join(wkspace,"rasterToPolygon\X_RasterToPolygon.shp")
field = "Value"

# Execute RasterToPolygon
arcpy.RasterToPolygon_conversion(inRaster, outPolygons, "NO_SIMPLIFY", field)


# Make a layer from the feature class
arcpy.MakeFeatureLayer_management(outPolygons,"RasterToPolygon_lyr")



# Process: Select Layer By Location
arcpy.SelectLayerByLocation_management(in_layer="RasterToPolygon_lyr", overlap_type="INTERSECT" ,\
                                       select_features=design_road, selection_type="NEW_SELECTION",invert_spatial_relationship= "NOT_INVERT")

# Write the selected features to a new featureclass
output_feature=os.path.join(arcpy.env.workspace,"X_selected_watersheds")
arcpy.CopyFeatures_management("RasterToPolygon_lyr",output_feature )



#Name: RandomPointsRandomValues.py
#Purpose: create random points with random values

outGDB = arcpy.env.workspace
outName = "X_random_points"
conFC = output_feature
numField = 200
arcpy.CreateRandomPoints_management(out_path=outGDB, out_name=outName, \
                                    constraining_feature_class=conFC, \
                                    number_of_points_or_field=numField, minimum_allowed_distance="1 meters")



# Name: FeatureToPolygon_Example2.py
# Description: Use FeatureToPolygon function to construct habitat areas
#              from park boundaries and rivers.

# Set local parameters
inFeatures = [output_feature, design_road]
outFeatureClass = os.path.join(arcpy.env.workspace,"X_featureToPolygon")
clusTol = "0.05 Meters"

# Use the FeatureToPolygon function to form new areas
arcpy.FeatureToPolygon_management(in_features=inFeatures,out_feature_class= outFeatureClass,
                                  label_features=os.path.join(outGDB,outName))




# Name: ZonalStatisticsAsTable_Ex_02.py
# Description: Summarizes values of a raster within the zones of 
#              another dataset and reports the results to a table.
# Requirements: Spatial Analyst Extension

# Import system modules
from arcpy.sa import *

# Set local variables
inZoneData = outFeatureClass
zoneField = "OBJECTID"
inValueRaster = SRTM_Digital_Elevation_Data_30m_tif
outTable = os.path.join(arcpy.env.workspace,"X_zonalStatsMean")


# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Execute ZonalStatisticsAsTable
outZSaT = ZonalStatisticsAsTable(inZoneData, zoneField, inValueRaster, 
                                 outTable, "NODATA", "MEAN")




# Name: AttributeSelection.py
# Purpose: Join a table to a featureclass and select the desired attributes

arcpy.env.qualifiedFieldNames = False
    
# Set local variables    
inFeatures = outFeatureClass # os.path.join(arcpy.env.workspace,"featureToPolygon")
layerName = "watershedMeanLyr" 
infield = "OBJECTID"
joinTable = outTable   # os.path.join(arcpy.env.workspace,"zonalStatsMean")
joinfield= "OBJECTID_1"
outFeature = os.path.join(arcpy.env.workspace,"X_watershed_with_mean")
    
# Create a feature layer from the vegtype featureclass
arcpy.MakeFeatureLayer_management (inFeatures,  layerName)
    
# Join the feature layer to a table
arcpy.AddJoin_management(in_layer_or_view=layerName, in_field=infield, join_table=joinTable,
                         join_field=joinfield, join_type="KEEP_COMMON" )
    
# Select desired features from veg_layer
#arcpy.SelectLayerByAttribute_management(layerName, "NEW_SELECTION", expression)
    
# Copy the layer to a new permanent feature class
arcpy.CopyFeatures_management(layerName, outFeature)


# Name: Project_Example.py
# Description: Projects a global coverage
# Requirements: ArcInfo Workstation



# projection 
# Local variables:
watershed_with_mean =  outFeature #os.path.join(arcpy.env.workspace,"X_watershed_with_mean")
watershed_with_mean_Project = os.path.join(arcpy.env.workspace,"watershed_with_mean_Project")

# Process: Project
arcpy.Project_management(watershed_with_mean, watershed_with_mean_Project, "PROJCS['WGS_1984_UTM_Zone_43N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',75.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", "", "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]", "NO_PRESERVE_SHAPE", "", "NO_VERTICAL")



# area calculation from watershed_with_mean_project feature class
import numpy as np

featureclass = watershed_with_mean_Project
field_names = ['OBJECTID', 'CID', 'MEAN', 'Shape_Area'] # [u'OBJECTID', u'Shape', u'CID', u'OBJECTID_1', u'OBJECTID_12', u'COUNT', u'AREA', u'MEAN', u'Shape_Length', u'Shape_Area']


# variables
da=[]
same=[]
area=[]
max_area=[]
area_selection_id=[]

with arcpy.da.SearchCursor(featureclass, field_names) as cursor:
    for row in cursor:
        da.append(row)
        same.append(row[1])  #CID

same=set(same)

aaa=np.array(da)

for i in same:
    area=[]
    for row in da:
        
        if row[1]==i and i!=0:# CID FIELD
            area.append(row)
    #print(i,area)
    if len(area)>1:
        stm=area[np.array(area)[:,2].argmax()] 
        max_area.append(stm[3]) # Appending area for each CID
        area_selection_id.append(stm[0]) # [OID]

summation=np.sum(max_area)
    
print('The area is: {} meter square'.format(summation))   



# selecting the contributing area

# local vairables
contributingArea = os.path.join(arcpy.env.workspace,"X_contributing_area")

# Make a layer from the feature class
arcpy.MakeFeatureLayer_management(watershed_with_mean_Project,"watershed_with_mean_ProjectLyr")

arcpy.SelectLayerByAttribute_management("watershed_with_mean_ProjectLyr", "NEW_SELECTION", ' "OBJECTID" in {} '.format(tuple(area_selection_id))) 

# Write the selected features to a new featureclass
arcpy.CopyFeatures_management("watershed_with_mean_ProjectLyr", contributingArea)


# cliping clustered area for calculating coef of runoff

#variables
#in_features=r"D:\6th_sem\final_year\arcGIS\data\clusters_rrzone.tif"
#clip_features = contributingArea  #os.path.join(arcpy.env.workspace,"contributing_area")
#out_feature_class = os.path.join(wkspace,"clustered_design_area\clustered_design_area_clip.tif")
#Use_Input_Features_for_Clipping_Geometry= "true"


# Process: Clip
# Set local variables
in_features = os.path.join(arcpy.env.workspace,"clusters_polygon_rrzone")
clip_features = contributingArea
X_clustered_design_area_clip = os.path.join(arcpy.env.workspace,"X_clustered_design_area_clip")
xy_tolerance = ""

arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("WGS 1984 UTM Zone 43N")
# Execute Clip
arcpy.Clip_analysis(in_features, clip_features, X_clustered_design_area_clip, xy_tolerance)


# raster to polygon













# Process: Create Random Points
#arcpy.CreateRandomPoints_management(currentWS_v2_gdb, Output_Point_Feature_Class, intersecting_parcels, "0 0 250 250", Number_of_Points__value_or_field_, Minimum_Allowed_Distance__value_or_field_, "POINT", "0")

# Process: Feature To Polygon
#arcpy.FeatureToPolygon_management("'karnataka_highway_Clip selection';'karnataka_highway_Clip selection'", parcels_divided_by_design_road, "", "ATTRIBUTES", random_points_on_parcels)

# Process: Zonal Statistics as Table
#arcpy.gp.ZonalStatisticsAsTable_sa(parcels_divided_by_design_road, Zone_field, SRTM_Digital_Elevation_Data_30m_tif, zonal_stats_elevation, "DATA", Statistics_type)

# Process: Join Field
#arcpy.JoinField_management(parcels_divided_by_design_road, Input_Join_Field, zonal_stats_elevation, Output_Join_Field, Join_Fields)




## Time of concentration ##

# taking the road length within contributing area for time of concentration
design_road_1_Project = design_road
contributing_area = contributingArea #os.path.join(arcpy.env.workspace,"X_contributing_area")
design_road_1_Project_Clip_C = os.path.join(arcpy.env.workspace,"X_design_road_1_Project_Clip_SaveLoc")  # gives the road network which intersects the contributng area
Points_on_Tc = os.path.join(arcpy.env.workspace,"X_PointsOnLineTc")
arcpy.GeneratePointsAlongLines_management(design_road_1_Project_Clip_C, Points_on_Tc, "DISTANCE", "100000 meters", "", "END_POINTS")

from arcpy.sa import *

# Set local variables
inZoneData = Points_on_Tc
zoneField = "OBJECTID"
inValueRaster = SRTM_Digital_Elevation_Data_30m_tif
outTable = os.path.join(arcpy.env.workspace,"X_ZonalStatsTc")


# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Execute ZonalStatisticsAsTable
X_zonalStats_Tc = ZonalStatisticsAsTable(inZoneData, zoneField, inValueRaster, 
                                 outTable, "NODATA", "MEAN")

featureclass1 = design_road_1_Project_Clip_C
featureclass2 = X_zonalStats_Tc
field_names1 = ['Shape_Length'] # for design_road_1_Project_Clip_C
field_names2 = ['MEAN'] # for X_zonalStats_Tc

# variables
tc_length=0
tc_mean=[]
with arcpy.da.SearchCursor(featureclass1, field_names1) as cursor_d:
    for row in cursor_d:
        tc_length=row[0]

with arcpy.da.SearchCursor(featureclass2, field_names2) as cursor1:
    for row in cursor1:
        tc_mean.append(row[0])

# slope
slope=(max(tc_mean)-min(tc_mean))/ (tc_length)
print('slope is: {}'.format(slope))
# sheet flow
# local variables
mannings_n= 0.011 # check
print('mannings_n is: {}'.format(mannings_n))
intensity_mm_hr = 9.0799  # mm
p2_24mm= intensity_mm_hr*24  # mm for 24 hr duration
slope_sheet= slope
length_sheet= 30 * 3.28 # 30 meters and converting to feet
Tc_sheetFlow = ((0.007*(mannings_n*length_sheet)**0.8)/((p2_24mm*0.0393701)**0.5 * (slope_sheet)**0.4))* 60 # in mins 

#channel flow
velocity_assumed= 2.0 # m/s
Tc_channel= (tc_length-30)/(60*velocity_assumed)
tc_total= Tc_channel+Tc_sheetFlow
print('tc_total: {}'.format(tc_total))
intensity_kotyari = (7.1*(2**0.2)*p2_24mm**0.33)/((tc_total/60)**0.71) # mm/hr

# coef of runoff
feature_class= X_clustered_design_area_clip
field_name= ['classType','Shape_Area']
coef_arearunoff=[] # calculate C1*A1 +C2*A2 ...
area_runoff= [] # calculate the total area A1+A2+...
classes_coef= {'Bareland':0.6 , 'Urban': 0.9,'Forest':0.2}  # check
with arcpy.da.SearchCursor(feature_class, field_name) as cursor:
    for row in cursor:
        coef_arearunoff.append(classes_coef[row[0]]*row[1])
        area_runoff.append(row[1])

coef_of_runoff= sum(coef_arearunoff)/sum(area_runoff)
print('coef_of_runoff is : {}'.format(coef_of_runoff))
# discharge
Q_dis= 0.028* coef_of_runoff * (intensity_kotyari/10) * sum(area_runoff)* 0.0001  # area in ha now, intensity in cm/hr
print('Discharge Q is: {}'.format(Q_dis))
# depth evaluation


# width
b_width= 1.5 # meters
depth_assumed= 2.0 # meters
side_slope = 0.0

# Python program for implementation 
# of Bisection Method forsolving equations 

def func(depth_assumed):

    A_f= (b_width +side_slope * depth_assumed)*depth_assumed  # in m2
    P_peri= b_width + 2*(side_slope**2+1)*depth_assumed  # in meters
    R_radius= float(A_f)/P_peri  #
    Q_cal= (1/mannings_n)*A_f*(R_radius**(0.666))*slope**(0.5)  # in m3/s
    return Q_dis-Q_cal

# Prints root of func(x) 
# with error of EPSILON 
def bisection(a,b): 

	if (func(a) * func(b) >= 0): 
		print("You have not assumed right a and b\n") 
		return

	c = a 
	while (abs(b-a) >= 0.01): 

		# Find middle point 
		c = (a+b)/2

		# Check if middle point is root 
		if (func(c) == 0.0): 
			break

		# Decide the side to repeat the steps 
		if (func(c)*func(a) < 0): 
			b = c 
		else: 
			a = c 
			
	print("The value of depth is : ","%.4f"%c) 
	
# Driver code 
# Initial values assumed 
a =20.0
b = 0.0
tc_total=bisection(a, b) 




























##design_road_1_Project = design_road
##contributing_area = contributingArea #os.path.join(arcpy.env.workspace,"X_contributing_area")
##design_road_1_Project_Clip_C = os.path.join(arcpy.env.workspace,"X_design_road_1_Project_Clip_SaveLoc")
##design_road_1_Project_Clip_G = r"D:\6th_sem\final_year\arcGIS\arcpy_currentWS.gdb\design_road_1_Project_Clip_G" #point  # point to be changed, calculated
##contributing_area__3_ = design_road_1_Project_Clip_G
##contributing_area__2_ = contributing_area
##contributing_area_PolygonToL =os.path.join(arcpy.env.workspace,"X_PolygonToLine_FromContributingArea") # "D:\\6th_sem\\final_year\\arcGIS\\arcpy_currentWS.gdb\\contributing_area_PolygonToL"
##contributing_area_PolygonToL1 = os.path.join(arcpy.env.workspace,"X_PointsOnLineFromContributingArea") #"D:\\6th_sem\\final_year\\arcGIS\\arcpy_currentWS.gdb\\contributing_area_PolygonToL1"
##X_pointDistanceTable = os.path.join(arcpy.env.workspace,"X_pointDistanceTable")  #"D:\\6th_sem\\final_year\\arcGIS\\arcpy_currentWS.gdb\\points_distance_table"

# Process: Clip (2)
#arcpy.Clip_analysis(design_road_1_Project, contributing_area, design_road_1_Project_Clip_C, "")

# make feature layer
#arcpy.MakeFeatureLayer_management(contributing_area__2_,"contributing_area__2_lyr")


# Process: Select Layer By Location
#arcpy.SelectLayerByLocation_management("contributing_area__2_lyr", "INTERSECT", design_road_1_Project_Clip_G, search_distance="1 meters, "NEW_SELECTION", "NOT_INVERT")

# Process: Polygon To Line
#arcpy.PolygonToLine_management("contributing_area__2_lyr", contributing_area_PolygonToL, "IGNORE_NEIGHBORS")

# Process: Generate Points Along Lines
#arcpy.GeneratePointsAlongLines_management(contributing_area_PolygonToL, contributing_area_PolygonToL1, "PERCENTAGE", "", "5", "END_POINTS")

# Process: Point Distance
#arcpy.PointDistance_analysis(design_road_1_Project_Clip_G, contributing_area_PolygonToL1, X_pointDistanceTable, "")


# joining the Points table with
                                       # actually no need to join
