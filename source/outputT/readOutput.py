## Obsolete script.
## Was used for debugging outputT / output to POVray.
## Use ~/__work/talks/standardPlots/tpmPlots/makeAnimation.py

## M.Mueller@sron.nl, 2014/05/30


#!/usr/bin/python
fileName = 'testoutput.dat'

file=open(fileName, 'r')
lines=file.readlines();
file.close()
# purge comment lines starting in '#', also eliminate trailing '\n'
lines=[line.strip() for line in lines if not line.startswith('#')]

# start parsing:
nTimes=int(lines[0])
shapeFileName = lines[1]
nFacets=int(lines[2])
# vectors to Sun:
dummy=lines[3].split(';')
assert len(dummy) == nTimes
A2SunX=[]
A2SunY=[]
A2SunZ=[]
for vector in dummy:
    bla=vector.split(',')
    assert len(bla)==3
    A2SunX.append(float(bla[0].strip()))
    A2SunY.append(float(bla[1].strip())) 
    A2SunZ.append(float(bla[2].strip()))
# vectors to observer:
dummy=lines[4].split(';')
assert len(dummy) == nTimes
A2ObsX=[]
A2ObsY=[]
A2ObsZ=[]
for vector in dummy:
    bla=vector.split(',')
    assert len(bla)==3
    A2ObsX.append(float(bla[0].strip()))
    A2ObsY.append(float(bla[1].strip())) 
    A2ObsZ.append(float(bla[2].strip()))
# JDs:
dummy=lines[5].split(';')
assert len(dummy) == nTimes
JDs=[]
for bla in dummy:
    JDs.append(float(bla.strip()))
# TSS
TSS=float(lines[6].strip())
# actual lightcurve in dimensionless temperature
uPerFacet=[]
for i in range(7, nFacets+7):
    bla=lines[i].split(';')
    assert i == int(bla[0].strip())+7
    bla=[float(x.strip()) for x in bla[1].split(',')]
    assert len(bla)==nTimes
    uPerFacet.append(bla)
assert len(uPerFacet)==nFacets
# Done parsing

# Read in shape model (strip '.convex' off of file used by TPM code)
import os, sys
shapesDir = os.path.join(os.path.expanduser('~'), '__work', 'shapes', 'spheres')
assert os.path.isdir(shapesDir)
sys.path.append(shapesDir)
from makePyramidNumpy import shapeModel
shape=shapeModel(shapeFileName[:shapeFileName.find('.convex')])
vertices=shape.vertices
facets=shape.facets
normalVectors=[]
# outward facing, norm ~ 1
for facet in facets:
    normal=vertices[facet[0]].copy()
    for i in [1,2]:
        normal += vertices[facet[i]]
    normalVectors.append(-normal/3.)
assert len(facets) == len(normalVectors)
assert len(facets) == nFacets

# output POVray files:
nColors=256 # nColor shades of blue/red
dummy=uPerFacet[0]
tMax=max(dummy)
for dummy in uPerFacet[1:]:
    if max(dummy) > tMax:
        tMax=max(dummy)
    pass
pass


## Coordinate system: shape models use right-hand system, with z pointing 'up.'
## POVray uses left-hand system, with y pointing 'up.'
## Solution: flip y and z in output.

for i in range(nTimes):
    temperatures=[lc[i] for lc in uPerFacet]
    assert len(temperatures) == nFacets
    out=open('step%i.pov'%i, 'w')
    out.write('#include "colors.inc"\n')
    out.write('camera {\n')
    out.write('       orthographic\n')
    #out.write('       location <%f,%f,%f>*20\n'%(A2ObsX[i], A2ObsY[i], A2Obzs[i]))
    out.write('       location <%f,%f,%f>*20\n'%(A2ObsX[i], A2ObsZ[i], A2ObsY[i])) #  flip y <--> z 
    out.write('       look_at <0,0,0>\n')
    out.write('       angle 20 // play with this until mesh fills image\n')
    out.write('}\n')
    out.write('\n')
    #out.write('light_source { <%f,%f,%f>*30 color White}\n'%(A2SunX[i], A2SunY[i], A2SunZ[i]))
    out.write('light_source { <%f,%f,%f>*30 color White}\n'%(A2SunX[i], A2SunZ[i], A2SunY[i])) #  flip y <--> z 
    out.write('\n')
    out.write('mesh2 { \n')
    # vertices
    out.write('\tvertex_vectors {\n')
    out.write('\t\t%i, // nVertices\n'%(len(vertices)))
    for vertex in vertices:
        out.write('\t\t<%f,%f,%f>, \n'%(vertex[0], vertex[2], vertex[1]))  #  flip y <--> z 
        #out.write('\t\t<%f,%f,%f>, \n'%(vertex[0], vertex[1], vertex[2]))
    out.write('\t}\n')
    # colors
    out.write('\ttexture_list {\n')
    out.write('\t\t%i, \n'%nColors)
    inc=1./float(nColors-1)
    for i in range(nColors):
        out.write('\t\ttexture {pigment{rgb<%f,0,%f>} finish {ambient 0.3} }\n'%(i*inc, (nColors-1-i)*inc))
        #out.write('\t\ttexture {pigment{rgb<%f,0,%f>} finish {specular 0.6 ambient 0.3} }\n'%(i*inc, (nColors-1-i)*inc))
        ## phong/specular was misleading: light appeared to come from wrong directions
    out.write('\t}\n')
    # facets
    out.write('\tface_indices {\n')
    out.write('\t\t%i, \n'%len(facets))
    for facet, t in zip(facets, temperatures):
        colorIndex = int(t/tMax*float(nColors-1))
        out.write('\t\t<%i,%i,%i>, %i, \n'%(facet[0], facet[1], facet[2], colorIndex))
    out.write('\t}\n')
    out.write('}\n')
    out.close()

#assert False

#nTimes=300

# run POVray on all of them:
# Run ~20 processes at once.  Not sure why it can't be more; OS appears to choke.
import subprocess, time
processes=[] # list of active processes
maxProcesses=20
nSleep=0
for i in range(nTimes):
    process=subprocess.Popen("nice povray -H320 -W400 step%i.pov"%i, shell=True)
    processes.append(process)
    while len(processes) == maxProcesses:
        # check if any processes are done, yet
        keepGoing=True
        for index, process in reversed(list(enumerate(processes))):
            # loop backwards so I can pop processes without messing up the indices
            if process.poll() is not None:
                processes.pop(index)
                keepGoing=False
        if keepGoing:
            # all processes still running; wait 0.1 sec before checking again
            time.sleep(0.1)
            nSleep += 1
    pass
# wait for all running processes to finish before moving on
for process in reversed(processes):
    if process.poll() is None:
        process.wait()
    processes.pop()
print "%i sleep cycles of 0.1 sec each"%nSleep

#assert False
cmd = "convert -delay 10 -loop 0 "
for i in range(nTimes):
    cmd+="step%i.png "%i
cmd += " animation.gif"
process=subprocess.Popen(cmd, shell=True)
process.wait()
