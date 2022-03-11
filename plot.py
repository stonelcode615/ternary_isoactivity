#! /usr/bin/python

from matplotlib import pyplot as plt
import sys
import optparse
import numpy as np


sin60=np.sin(np.pi/3)

def file_reading(filename):
    file_handle = open(filename,'r')
    index = 0
    row_item = []
    while True:
        s1 = file_handle.readline()
        if s1 == '':
            break
        else:
            s2 = s1.rstrip('\n').split()
        try:
            float(s2[0])
        except ValueError:
            topname = s2[:]
            continue
        content = [float(item) for item in s2]
        row_item.append(content)
    return topname, row_item


def transform1(RawData):
    DisplayData=[]
    for coordinate in RawData:
        NewData=[]
        NewData.append(coordinate[0]+0.5*coordinate[1])
        NewData.append(coordinate[1]*sin60)
        DisplayData.append(NewData)
    return DisplayData

def transform2(RawData):
    DisplayData=[]
    for coordinate in RawData:
        NewData=[]
        NewData.append(coordinate[0]+0.5*coordinate[1])
        NewData.append(coordinate[1]*sin60)
        NewData.append(coordinate[2]+0.5*coordinate[3])
        NewData.append(coordinate[3]*sin60)
        DisplayData.append(NewData)
    return DisplayData

def binodal_plot(ComponentName,coordinate,frame,pn):
    symbols=['b-','g--']
    legend=['Binodal','Spinodal']

    x_alpha = []
    x_beta  = []
    y_alpha = []
    y_beta  = []
    for num,item in enumerate(coordinate):
        x_alpha.append(item[0])
        y_alpha.append(item[1])
        x_beta.append(item[2])
        y_beta.append(item[3])
    frame.plot(x_alpha,y_alpha,symbols[pn],label=legend[pn])
    frame.plot(x_beta,y_beta,symbols[pn])


def exp_tieline_plot(ComponentName,coordinate,frame):
    for num, item in enumerate(coordinate):
        if num==0:
            frame.plot([item[0],item[2]],[item[1],item[3]],'ro-',label='Exp')
        else:
            frame.plot([item[0],item[2]],[item[1],item[3]],'ro-')
    


def fit_tieline_plot(ComponentName,coordinate,frame,pn):
    symbols=['ks:','bs--','g^:','k>:']
    legend='Fit %d' %pn
    for num, item in enumerate(coordinate):
        if pn==0:
            if num==0:
                frame.plot([item[0],item[2]],[item[1],item[3]],symbols[pn],mfc='None',label='Fit')
            else:
                frame.plot([item[0],item[2]],[item[1],item[3]],symbols[pn],mfc='None')
        else:
            if num==0:
                frame.plot([item[0],item[2]],[item[1],item[3]],symbols[pn],label=legend)
            else:
                frame.plot([item[0],item[2]],[item[1],item[3]],symbols[pn])
    
def ticks_plot(vertices, tick_number,frame):
    tick_scale = 0.2
    ticks_vector=[]
    ticks_vector.append(tick_scale*(vertices[0]-vertices[2])/tick_number)
    ticks_vector.append(tick_scale*(vertices[1]-vertices[0])/tick_number)
    ticks_vector.append(tick_scale*(vertices[2]-vertices[1])/tick_number)
    segments = np.linspace(0,1,tick_number+1)
    for i in range(3):
        if i+1==3:
            j = 0
        else:
            j = i+1
        x = (vertices[j][0]-vertices[i][0])*segments+vertices[i][0]
        x = np.vstack((x,x+ticks_vector[i][0]))
        y = (vertices[j][1]-vertices[i][1])*segments+vertices[i][1]
        y = np.vstack((y,y+ticks_vector[i][1]))
        plt.plot(x,y,'k',lw=1)
        label_fraction = [str(index/5.) for index in range(1,5)]
        for k in range(1,5):
            if i==0:
                frame.text(x[1][k*2]-0.03,y[1][k*2]-0.04, label_fraction[k-1])
            elif i==1:
                frame.text(x[1][k*2]+0.017,y[1][k*2]-0.01,label_fraction[k-1])
            else:
                frame.text(x[1][k*2]-0.06,y[1][k*2]+0.01,label_fraction[k-1])

def grids_plot(vertices, grid_number, frame):
    segments = np.linspace(0,1,grid_number+1)
    x_coordinate=[]
    y_coordinate=[]
    for i in range(3):
        j = i+1
        if j==3:
            j=0
        x_coordinate.append((vertices[j][0]-vertices[i][0])*segments+vertices[i][0])
        y_coordinate.append((vertices[j][1]-vertices[i][1])*segments+vertices[i][1])
    for i in range(3):
        j = i+1
        if j==3:
            j=0
        x = np.vstack((x_coordinate[i],x_coordinate[j][::-1]))
        y = np.vstack((y_coordinate[i],y_coordinate[j][::-1]))
        plt.plot(x,y,'k:',lw=1)



def triangle_plot(points,sides,ticks_on, grids_on, ticks_number, frame):
    frame.plot(sides[0,1:],sides[1,1:],'k-',lw=2.0)
    frame.plot(sides[0,:2],sides[1,:2],'k-',lw=3.75)
    if grids_on == 'on':
        grids_plot(points, ticks_number, frame)
    if ticks_on == 'on':
        ticks_plot(points, ticks_number, frame)
    


def make_format(current,composition):
    def format_coord(x,y):
        display_coord = current.transData.transform((x,y))
        inv = current.transData.inverted()
        ax_coord = inv.transform(display_coord)
        coords = ax_coord
        xy = []
        xy.append(coords[0]-coords[1]/np.sqrt(3.0))
        xy.append(coords[1]/sin60)
        xy.insert(0,1-xy[0]-xy[1])
        string0 = '{0[0]}: {1[0]:.2%}   {0[1]}: {1[1]:.2%}  {0[2]}: {1[2]:.2%}'.format(composition,xy)
        print string0
        return string0
    return format_coord

#################################################

vertex1 = np.array([0,0])
vertex2 = np.array([1,0])
vertex3 = np.array([0.5,sin60])
vertices = np.vstack((vertex1, vertex2, vertex3))
edges = np.concatenate((vertex1.reshape(2,1),vertex2.reshape(2,1),vertex3.reshape(2,1),vertex1.reshape(2,1)),axis=1)

plt.ion()
fig = plt.figure(figsize=(8,6),dpi=80,facecolor='w')
ax = plt.subplot(111)
ax.axis('equal')
ax.axis('off')

parser = optparse.OptionParser()
parser.add_option('--tcal',dest='Cal_tieline',help='Binodal Plot',action="append")
parser.add_option('--texp',dest='Exp_tieline',help='TieLine Plot')
parser.add_option('--tfit',dest='Fit_tieline',help='Tie Test Plot',action="append")

parser.add_option('--grid',dest='GridsOn',help='plot grids, on/off')
parser.add_option('--tick',dest='TicksOn',help='plot ticks, on/off')

(options, args) = parser.parse_args()

triangle_plot(vertices,edges,options.TicksOn,options.GridsOn,10,ax)


component_name=[]
if options.Cal_tieline is None:
    print     'No calculated tielines to Plot'
else:
    for fn,fname in enumerate(options.Cal_tieline):
        component, row = file_reading(fname)
        component_name = component
        display_coordinate = transform2(row)
        binodal_plot(component, display_coordinate,ax,fn)

if options.Exp_tieline is None:
    print     'No experital  tielines to Plot'
else:
    component, row = file_reading(options.Exp_tieline)
    component_name = component
    display_coordinate = transform2(row)
    exp_tieline_plot(component, display_coordinate,ax)

if options.Fit_tieline is None:
        print 'No fitting    tielines to Plot'
else:
    for fn,fname in enumerate(options.Fit_tieline):
        component, row = file_reading(fname)
        component_name = component
        display_coordinate = transform2(row)
        fit_tieline_plot(component, display_coordinate,ax,fn)



ax.text(-0.05,-0.05,component[0],fontweight='bold')
ax.text(1.05,-0.05, component[1],fontweight='bold')
ax.text(0.42,0.90,  component[2],fontweight='bold')

ax.format_coord = make_format(ax,component_name)
print component
#ax.format_coord = lambda x,y: '%s:%.2f   %s:%.2f   %s:%.2f'%(component[0],1-y/sin60-(x-y)/np.sqrt(3.0),component[1],(x-y)/np.sqrt(3.0),component[2],y/sin60)
ax.format_coord = lambda x,y: '%s:%.2f   %s:%.2f   %s:%.2f'%(component[0],round(np.absolute(x+np.sqrt(3.0)/3.*y-1.0),2),component[1],round(np.absolute(x-np.sqrt(3.0)/3.*y),2),component[2],round(2*y/np.sqrt(3.0),2))

plt.legend()
plt.show()

plt.savefig('ternary.eps')
exit=raw_input('enter any key to exit')



