# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 14:06:13 2019
Exports reach polygon shapefiles given input of centerlines
add watershed boundaries:
@author: john
"""

# Import the Earth Engine Python Package
import time
import ee

# Initialize the Earth Engine object, using the authentication credentials.
ee.Initialize()

def StringToNum(feature):
    num = ee.Number.parse(feature.get('reach_id')).toInt64()
    return feature.set('reach_ID', num)

### define functions for finding pixels  closest to each reach in centerline file
# change the "ID" to whatever the unique ID is of each reach
def fdtFun(f):
    dt = ee.Image.constant(0).byte().paint(ee.FeatureCollection(f), 1).fastDistanceTransform()
    return dt.copyProperties(f, ['reach_ID'])

def mfun(c, p):
    return ee.Image(c).min(p)

def minFun(f):
  return ee.Image(f.iterate(mfun, ee.Image(f.first())))

def idFun(c, p):
    c = ee.Image(c)
    p = ee.Image(p)
    return p.where(c.eq(minMap), ee.Number(c.get('reach_ID')))

# functions for bacth processing
def maximum_no_of_tasks(MaxNActive, waitingPeriod):
  ##maintain a maximum number of active tasks
  time.sleep(10)
  ## initialize submitting jobs
  ts = list(ee.batch.Task.list())

  NActive = 0
  for task in ts:
       if ('RUNNING' in str(task) or 'READY' in str(task)):
           NActive += 1
  ## wait if the number of current active tasks reach the maximum number
  ## defined in MaxNActive
  while (NActive >= MaxNActive):
      time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait for 2min and check again
      ts = list(ee.batch.Task.list())
      NActive = 0
      for task in ts:
        if ('RUNNING' in str(task) or 'READY' in str(task)):
          NActive += 1
  return()
  
def getCols(tableMetadata):
  print(tableMetadata.columns)


#load centerlines and convert reach from string to number
lines = ee.FeatureCollection("users/johngardner87/SWORD_SA").map(StringToNum)

# load amazon basin
hb3 = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_3")

hb3_amazon = hb3.filter(ee.Filter.eq('PFAF_ID', 622))

#print(type(hb3_amazon))

lines_clean = (lines.filterBounds(hb3_amazon.geometry()))\
.filter(ee.Filter.lt('type', 4))

print(lines_clean.limit(1).getInfo())

hb6 = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_6")

hb6_amazon = hb6.filterBounds(hb3_amazon.geometry())

hbsort = hb6_amazon.sort('PFAF_ID')

hbID = hbsort.aggregate_array('PFAF_ID').getInfo()

# make a folder in your google drive manually to output data
#dlDir = 'D:/GoogleDrive/SWORD_polygons_amazon' 
#filesDown = os.listdir(dlDir)  # -->
#filesDown = [int(i.replace(".shp", "")) for i in filesDown]
#hbID  = [i for i in hbIF if i not in filesDown]

#for x in range(0,len(hbID)):
    
for x in range(0,2):
    
    #hb = hb6_amazon.filterMetadata('PFAF_ID' , 'equals', hbID[x])
    hb = hb6_amazon.filter(ee.Filter.eq('PFAF_ID' ,  hbID[x]))
    
    #print(hb.limit(1).getInfo())
    
    ### THIS LINE FAILS
    linesClip = lines_clean.filterBounds(hb.geometry())
    #linesClip = lines_clean.filter(ee.Filter.intersects('.geo', hb))
    
    print(linesClip.limit(1).getInfo())

# set maximum buffer from center line to assign pixels a reach ID (10,000 m)
    linesClip_buffer = linesClip.geometry().buffer(10000)

    #print(linesClip_buffer.getInfo())
    
    fdt = linesClip.map(fdtFun)
    
    minMap = minFun(fdt)
    
    idImage = fdt.iterate(idFun, ee.Image(fdt.first()))
    
    print(idImage.getInfo())
### make image and clip to max distance buffer from centerline for vector export
    #reachID = ee.Image(idImage).rename('reach_id').clip(linesClip_buffer)
    reachID = ee.Image(idImage).clip(linesClip_buffer)
    
    print(reachID.getInfo())
##
    vectors = reachID.addBands(reachID).reduceToVectors(
        crs = reachID.projection(), 
        scale =  90, 
        geometry = hb, 
        geometryType = 'polygon',
        eightConnected= False,
        tileScale = 16,
        bestEffort = True,
        labelProperty = 'reach_ID',
        maxPixels = 3000000000000000, 
        reducer = ee.Reducer.first()
        )

    
    dataOut = ee.batch.Export.table.toDrive(collection = vectors, \
                                              description = str(hbID[x]),\
                                              folder = 'SWORD_polygons_amazon',\
                                              fileFormat = 'shp')
    
    print(hbID[x])
    
    
  #Check how many existing tasks are running and take a break if it's >15  
    maximum_no_of_tasks(15, 60)
  #Send next task.
    dataOut.start()

#Make sure all Earth engine tasks are completed prior to moving on.  
maximum_no_of_tasks(1,300)
print('done')

### NOTES
## watersheds missed due to exceeeded memory = 1701, 1704, 1706, 1018, 1101, 0708, 0601

### Scratch text
#    fdt = nhdClip.map(function(f) {
#    dt = ee.Image.constant(0).byte().paint(ee.FeatureCollection(f), 1).fastDistanceTransform()
#    return(dt.copyProperties(f, ['COMID']))})

#    minMap = ee.Image(fdt.iterate(function(c, p) {return(ee.Image(c).min(p))}, ee.Image(fdt.first())))

#    idImage = fdt.iterate(function(c, p) {
#    c = ee.Image(c)
#    p = ee.Image(p)
#    return(p.where(c.eq(minMap), ee.Number(c.get('COMID'))))},
#    ee.Image(fdt.first()))
