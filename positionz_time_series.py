#!/bin/python

import sys
import numpy as np
import os.path
import datetime
import re
from collections import namedtuple
from ellipsoid import grs80
from functools import partial

tsfields='name epoch x y z'.split()
tsfilename='{code}_igs08_xyz.dat'

base_dir = os.path.dirname(__file__)
stn_file='stations.txt'
obs_file='h1_coords.txt'

tsobs=namedtuple('tsobs','code time xyz')
xyzfile=namedtuple('xyzfile','code filename')

def robust_standard_error(obs):
    errors=[0]*3
    for axis in range(3):
        diffs=np.abs(obs[1:,axis]-obs[:-1,axis])
        se=np.percentile(diffs,95.0)/(1.96*np.sqrt(2))
        errors[axis]=se
    return errors

def readdata(filename,code=None,xyz0=None,enu=False,sort=False,transform=None,enu_factor=1.0):
    '''
    Reads timeseries data from the specified file.  The file is assumed to contain whitespace
    delimited data.  Each line contains the station code, the date/time in format YYYY-MM-DDThh:mm:ss
    where T is a literal character, and the X, Y, and Z coordinates, and an optional text data field.  
    The first line of the file must be a header containing 'name epoch x y z'.
    
    The function call can provide the following optional parameters:

        code    Defines the code expected - data not matching are rejected
        xyz0    Defines a reference coordinate for ENU coordinates - default is the first coordinate
        enu     If true then data are returned as ENU relative to reference rather than XYZ
        sort    If true then data are sorted by date, default is the order in the file
        transform   An optional function to apply a transformation function to the XYZ data.  This
                    is called as xyz=transform(xyz,date)
        enu_factor  A scale factor to apply to ENU data
    '''
    if not os.path.exists(filename):
        raise RuntimeError('Timeseries observation file '+of+' doesn\'t exist')
    obs=[]
    first=True
    rmat=None
    with open(filename) as f:
        fields=f.readline().lower().split()
        if len(fields) < len(tsfields) or fields[:len(tsfields)] != tsfields:
            raise RuntimeError('Station file '+of+' doesn\'t have the correct fields')
        for l in f:
            parts=l.split()
            if len(parts) < 5:
                continue
            fcode=parts[0].upper()
            xyz = np.array([float(p) for p in parts[2:5]])
            otime = datetime.datetime.strptime(parts[1],'%Y-%m-%dT%H:%M:%S')
            if code and fcode != code:
                continue
            if transform:
                xyz=transform(xyz,otime)
            if first:
                first=False
                if not xyz0:
                    xyz0=xyz
                if enu:
                    lon,lat,hgt=grs80.geodetic(xyz)
                    rmat=grs80.enu_axes(lon,lat)
            if enu:
                xyz=[(x-x0)*enu_factor for x,x0 in zip(xyz,xyz0)]
                xyz=rmat.dot(xyz)
            o=tsobs(fcode,otime,xyz)
            if sort:
                obs.append(o)
            else:
                yield tsobs(fcode,otime,xyz)
    if sort:
        obs.sort(key=lambda o: o.time)
        for o in obs:
            yield o

def finddata(directory,filename=tsfilename,sortby='code'):
    '''
    Retrieves a list of timeseries files from the specified directory.

    Returns a list of tsdata objects.  Takes an optional filename parameter
    which is the expected name of the time series data file including 
    {code} where the four character station code is expected.  Also takes
    a sortby argument which can be 'code' or 'z'.
    '''
    if '{code}' not in filename:
        raise RuntimeError('XYZ filename must contain {code} in xyzfiles')
    if not os.path.isdir(directory):
        raise RuntimeError(directory+' is not a directory in xyzfiles')
    f0,f1=filename.split('{code}',1)
    fre=re.compile('^'+re.escape(f0)+r'(\w{4})'+re.escape(f1)+'$')
    files=[]
    for fn in os.listdir(directory):
        m=fre.match(fn)
        if not m:
            continue
        code=m.group(1).upper()
        filename=os.path.join(directory,fn)
        try:
            files.append(tsdata(filename,code))
        except:
            raise
            pass
    if str(sortby).lower() == 'code':
        files.sort(key=lambda f: f.code)
    elif str(sortby).lower() == 'z':
        files.sort(reverse=True,key=lambda f: f.xyz[2])
    return files

class tsdata:

    def __init__(self,filename,code=None):
        data=None
        for l in readdata(filename,code):
            data=l
            break
        if not data:
            raise RuntimeError("No timeseries data in file "+filename)
        self.filename=filename
        self.code = data.code
        self.xyz=data.xyz
        self.llh=grs80.geodetic(data.xyz)

    def data( self, **options ):
        for d in readdata(self.filename,**options):
            yield d

    def data_array(self,enu=False,factor=0):
        '''
        Return data as an array of dates and an array of xyz coordinates, ie

        dates,xyzdata=f.data_array()
        '''
        dates=[]
        obs=[]
        for o in self.get_obs(enu,factor,True):
            dates.append(o.time)
            obs.append(o.xyz)
        return dates,np.array(obs)

