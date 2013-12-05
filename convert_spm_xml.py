import sys
import os.path
import os
import re
import argparse
from datetime import datetime, timedelta
from datetime import date as datetype
from collections import namedtuple
import numpy as np
from scipy.special import erf
from scipy.optimize import leastsq
import math
import ellipsoid
from xml.etree import ElementTree
from xml.dom import minidom


refdate=datetime(2000,1,1)
daysperyear=365.25

# Time conversion/utility functions

def asday( date ):
    if type(date) == float:
        return date
    if type(date)==str:
        try:
            date=datetime.strptime(date,"%d-%m-%Y %H:%M")
        except:
            try:
                date=datetime.strptime(date,"%d-%m-%Y")
            except:
                date=datetime.strptime(date,'%Y-%m-%dT%H:%M:%S')
    td = date.replace(tzinfo=None) - refdate
    return float(td.days)+float(td.seconds)/(60*60*24) 

def fromday( days ):
    if type(days) == datetime:
        return days
    td = timedelta(days)
    return refdate+td

def days_array( dates ):
    if not isinstance(dates,np.ndarray):
        if not isinstance(dates,list):
            dates=[dates]
        dates=np.array(dates)
    if type(dates[0]) in (datetime, datetype):
        dates=np.array([asday(d) for d in dates])
    return dates


# Parameter types used to define model parameters

class parameter( object ):

    def __init__(self,model,code,name,index,factor=1.0,format='{0}'):
        self._model=model
        self._code=code
        self._name=name
        self._index=index
        self._factor=factor
        self._format=format
        self._error=0.0
        self._calcdate=None
        self._isLinear=False
        self._covarIndex=-1
        self._saved=None

    def code( self ):
        return self._code

    def name( self ):
        return self._name

    def fixed( self ):
        return self._model._fixed[self._index]

    def setFixed( self, fixed ):
        self._model._fixed[self._index]=bool(fixed)

    def setValue( self, valuestr ):
        value=float(valuestr)
        value /= self._factor
        self.setFitValue(value)
        self._covarIndex=-1
        self._calcdate=None

    def getValue( self ):
        value=self.fitValue()
        value *= self._factor
        return self._format.format(value)

    def getError( self ):
        if self._error is None:
            return ''
        return self._format.format(self._error*self._factor)

    def covarIndex( self ):
        return self._covarIndex

    def calcDate( self ):
        return self._calcdate

    def toXmlValue( self ):
        return self.getValue()

    def fromXmlValue( self, xmlstr ):
        self.setValue( xmlstr )

    def fitValue( self ):
        return self._model._param[self._index]

    def setFitValue( self, value, index=-1, error=None ):
        self._model._param[self._index] = value
        # Assume that if error is provided it has been calculated, so reset calc date
        self._error=error
        if index >= 0 and error != None:
            self._covarIndex=index
            self._calcdate=datetime.now()

    def saveValue( self ):
        self._saved=[self.fitValue(),self._error,self._calcdate,self._covarIndex]

    def restoreValue( self ):
        if self._saved is not None:
            value,self._error,self._calcdate,self._covarIndex = self._saved
            self.setFitValue(value)

    def xmlElement( self ):
        element=ElementTree.Element('parameter')
        element.set('code',self.code())
        element.set('value',self.toXmlValue())
        element.set('fit','no' if self.fixed() else 'yes')
        if self._error:
            element.set('error',self.getError())
        if self._calcdate != None:
            element.set('calc_date',self._calcdate.strftime('%Y-%m-%dT%H:%M:%S'))
        if self._covarIndex >= 0:
            element.set('covariance_index',str(self._covarIndex))
        return element

    def loadXmlElement( self, element ):
        assert(element.tag=='parameter')
        assert(element.get('code')==self.code())
        self.fromXmlValue(element.get('value',self.toXmlValue()))
        self.setFixed(element.get('fit','yes') != 'yes')
        self._error=None
        self._calcdate=None
        self._covarIndex=-1
        error=element.get('error','')
        calcdate=element.get('calc_date','')
        covarIndex=element.get('covariance_index','')
        if error:
            self._error = float(error)/self._factor
        if calcdate:
            try:
                self._calcdate=datetime.strptime(calcdate,'%Y-%m-%dT%H:%M:%S')
            except:
                pass
        if covarIndex:
            try:
                self._covarIndex=int(covarIndex)
            except:
                pass

    def __str__( self ):
        return self.getValue()

class offset_parameter( parameter ):
    def __init__(self,model,code,name,index):
        parameter.__init__( self,model,code,name,index,factor=1000,format="{0:.1f}")
        self._isLinear=True
        
class velocity_parameter( parameter ):
    def __init__(self,model,code,name,index):
        parameter.__init__( self,model,code,name,index,factor=1000*365.25,format="{0:.3f}")
        self._isLinear=True


class date_parameter( parameter ):
    def __init__(self,model,code,name,index,format="{0:.1f}"):
        parameter.__init__(self,model,code,name,index)

    def setValue(self, valuestr ):
        self.setFitValue(asday(valuestr))

    def getValue(self):
        return fromday(self.fitValue()).strftime("%d-%m-%Y %H:%M")

    def getError( self ):
        return '' if self._error is None else "{0:.2f}".format(self._error)

    def toXmlValue(self):
        return fromday(self.fitValue()).strftime('%Y-%m-%dT%H:%M:%S')


# Functions used to build model

class base_function( object ):

    def __init__( self, model, date, nparam, params=None ):
        if not date:
            date = model.refdate
        self.model = model
        self.parameters=[]
        self._event=''
        self._nparam = nparam+1
        self._fixed = [False] * (nparam+1)
        self._param = [0.0] * (nparam+1)
        self._param[nparam] = asday(date)
        self._fixed[nparam] = True
        self._funcFixed = True
        self._enabled = True
        if params:
            self._param[:nparam] = params[:nparam]

    def setModel( self, model ):
        self.model = model

    def eventDate( self ):
        return fromday(self._param[-1])

    def _dateOffset( self, dates ):
        d = days_array(dates)-self._param[-1]
        return d.reshape(d.size,1)

    def setComponent( self, i, value, fixed=False ):
        assert i >= 0 and i < self._nparam-1
        self._param[i] = float(value)
        self._fixed[i] = bool(fixed)

    def fixed(self):
        return self._funcFixed

    def setFixed( self, fixed):
        self._funcFixed = fixed

    def enabled( self ):
        return self._enabled

    def setEnabled( self, enabled ):
        self._enabled = enabled

    def setEventName( self, name, detail='' ):
        self._event=name

    def eventName( self ):
        return self._event

    def eventDetail( self ):
        return self._eventDetail

    def xmlElement( self ):
        element=ElementTree.Element('component')
        element.set('type',type(self).__name__)
        if self._event:
            element.set('name',self._event)
        if not self._enabled:
            element.set('exclude','yes')
        element.set('fit','no' if self._funcFixed else 'yes')
        for p in self.parameters:
            element.append(p.xmlElement())
        return element

    @staticmethod
    def fromXmlElement( model, element ):
        classname = element.tag
        def getsubclass( base, classname ):
            for c in base.__subclasses__():
                if c.__name__==classname:
                    return c
                subc=getsubclass(c,classname)
                if subc:
                    return subc
        cclass=getsubclass(base_function,classname)
        if not cclass:
            raise ValueError('Invalid station model type '+classname)
        component=cclass(model)
        component._event=element.get('name','')
        component._enabled=element.get('exclude','No').lower() != 'yes'
        component._funcFixed=element.get('fit','No').lower() != 'yes'
        for pelement in element:
            code=pelement.get('code','')
            for p in component.parameters:
                if p.code() == code:
                    p.loadXmlElement(pelement)
                    break
        return component

    def paramStr( self ):
        return (self.eventDate().strftime('%d-%m-%Y') + ' [' +
            ', '.join([str(p) for p in self._param[:-1]]) + ']')

    def __str__( self ):
        strval = self._event + ': ' if self._event else ''
        return strval + self.paramStr()
    
class offset( base_function ):

    def __init__( self, model, offset=None ):
        base_function.__init__( self, model, refdate, 3, offset )
        self.parameters=[
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]
        self._funcFixed = False

    def calc( self, date ):
        return np.array([self._param[:3]]*len(date))

    def paramStr( self ):
        return type(self).__name__+ '  [{0:.1f}, {1:.1f}, {2:.1f}] mm'.format(
            *[x*1000 for x in self._param[:3]])

class velocity( base_function ):

    def __init__( self, model, velocity=None ):
        base_function.__init__( self, model, None, 3, velocity )
        self.parameters=[
            velocity_parameter(self,'ve_mmpy','East (mm/yr)',0),
            velocity_parameter(self,'vn_mmpy','North (mm/yr)',1),
            velocity_parameter(self,'vu_mmpy','Up (mm/yr)',2),
            ]
        self._funcFixed = False

    def calc( self, date ):
        y=self._dateOffset(date)
        return y.dot([self._param[:3]])

    def paramStr( self ):
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm/year'.format(
            *[x*daysperyear*1000 for x in self._param[:3]])

class velocity_change( base_function ):

    def __init__( self, model, date=refdate, velocity_change=None ):
        base_function.__init__( self, model, date, 3, velocity_change )
        self.parameters=[
            date_parameter(self,'date','Date of change',3),
            velocity_parameter(self,'ve_mmpy','East change (mm/yr)',0),
            velocity_parameter(self,'vn_mmpy','North change (mm/yr)',1),
            velocity_parameter(self,'vu_mmpy','Up change (mm/yr)',2),
            ]

    def calc( self, date ):
        y=np.maximum(self._dateOffset(date),0.0)
        return y.dot([self._param[:3]])

    def paramStr( self ):
        f=daysperyear*1000
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm/year at {3:%d-%m-%Y}'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, self.eventDate())

class cyclic( base_function ):

    def __init__( self, model, frequency, cosp=None, sinp=None ):
        if cosp and sinp:
            params = [cosp[0],cosp[1],cosp[2],sinp[0],sinp[1],sinp[2]]
        else:
            params=None
        base_function.__init__( self, model, refdate, 6, params )
        self._frequency=frequency
        self.parameters=[
            offset_parameter(self,'ecos_mm','East cosine (mm)',0),
            offset_parameter(self,'esin_mm','East sine (mm)',3),
            offset_parameter(self,'ncos_mm','North cosine (mm)',1),
            offset_parameter(self,'nsin_mm','North sine (mm)',4),
            offset_parameter(self,'ucos_mm','Up cosine (mm)',2),
            offset_parameter(self,'usin_mm','Up sine (mm)',5),
            ]

    def calc( self, date ):
        y=(self._dateOffset(date)/daysperyear)*math.pi*2*self._frequency;
        return np.cos(y).dot([self._param[:3]])+np.sin(y).dot([self._param[3:6]])

    def setComponent( self, i, issine, value, fixed=False ):
        assert type(issine) == bool, repr(issine)
        assert i >= 0 and i < 3
        if issine:
            i = i+3
        base_function.setComponent(self,i,value,fixed)

    def paramStr( self ):
        f=daysperyear*1000
        return (type(self).__name__+ 
                '  [{0:.2f}*cos+{3:.2f}*sin, {1:.2f}*cos+{4:.2f}*sin, {2:.2f}*cos+{5:.2f}*sin] mm'.format(
                    *[x*1000 for x in self._param[:6]])+
                ' frequency {0:.0f}/year'.format(self._frequency)
               )

class annual( cyclic ):

    def __init__(self,model, cosp=None,sinp=None):
        cyclic.__init__(self,model, 1.0,cosp,sinp)

class semiannual( cyclic ):

    def __init__(self,model, cosp=None,sinp=None):
        cyclic.__init__(self,model, 2.0,cosp,sinp)

class equipment_offset( base_function ):

    def __init__( self, model, date=refdate, offset=None ):
        base_function.__init__( self, model, date, 3, offset )
        self.parameters=[
            date_parameter(self,'date','Date of offset',3),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def calc( self, date ):
        y=np.where(self._dateOffset(date)>0,1,0)
        return y.dot([self._param[:3]])

    def paramStr( self ):
        f=1000;
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm at {3:%d-%m-%Y}'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, self.eventDate())

class tectonic_offset( equipment_offset ):

    def __init__( self, model, date=refdate, offset=None ):
        equipment_offset.__init__( self, model, date, offset )
        self.parameters[0]=date_parameter(self,'date','Date of event',3)

class slow_slip( base_function ):

    def __init__( self, model, date=refdate, duration=10.0/3.92, offset=None ):
        base_function.__init__( self, model, date, 4 )
        self.setDuration(duration)
        if offset:
            self._param[:3]=offset[:3]
        self.parameters=[
            date_parameter(self,'mid_date','Central date of event',4),
            parameter(self,'duration_days','Duration of event (95% of slip) (days)',3,factor=3.92,format='{0:.1f}'),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def setDuration( self, duration, fixed=True ):
        self.setComponent(3,duration,fixed)

    def calc( self, date ):
        y=0.5*(1+erf(self._dateOffset(date)/self._param[3]))
        # y=np.minimum(1,np.maximum(0,self._dateOffset(date)/self._param[3]))
        return y.dot([self._param[:3]])
    
    def paramStr( self ):
        f=1000;
        start=self.eventDate()
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm centred {3:%d-%m-%Y} 95% width {4:.2f} days'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, start, 3.92*self._param[3])

class slow_slip_ramp( base_function ):

    def __init__( self, model, date=refdate, end_date=None, offset=None ):
        base_function.__init__( self, model, date, 4 )
        if end_date is None:
            end_date=self._param[4]+10.0
        self.setEndDate(end_date)
        if offset:
            self._param[:3]=offset[:3]
        self.parameters=[
            date_parameter(self,'start_date','Start date of event',4),
            date_parameter(self,'end_date','End date of event',3),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def setEndDate( self, end_date, fixed=True ):
        self.setComponent(3,asday(end_date),fixed)

    def calc( self, date ):
        duration=self._param[3]-self._param[4]
        y=self._dateOffset(date)/duration
        y=np.minimum(1.0,np.maximum(y,0.0))
        return y.dot([self._param[:3]])

    def paramStr( self ):
        f=1000;
        start=self.eventDate()
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm between {3:%d-%m-%Y} and {4:%d-%m-%Y}'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, start, fromday(self._param[3]))

class exponential_decay( base_function ):

    def __init__( self, model, date=refdate, decay=10.0/math.log(2.0), offset=None ):
        base_function.__init__( self, model, date, 4 )
        self.setDuration(decay)
        if offset:
            self._param[:3]=offset[:3]
        self.parameters=[
            date_parameter(self,'date','Start date of event',4),
            parameter(self,'halflife_days','Half life of event (days)',3,factor=math.log(2.0),format='{0:.1f}'),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def calc( self, date ):
        y=self._dateOffset(date)/self._param[3]
        y=np.maximum(0,1.0-np.exp(-y))
        return y.dot([self._param[:3]])

    def setDuration( self, decay, fixed=True ):
        self.setComponent(3,decay,fixed)

    def paramStr( self ):
        f=1000;
        start=self.eventDate()
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm exponential decay from {3:%d-%m-%Y} half life {4:.1f} days'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, start, self._param[3]*math.log(2.0))

class exclude_obs( object ):

    def __init__( self, date, comment="", index=-1 ):
        self.date=date
        self.day=asday(date)
        self.index=index
        self.comment=comment

    def xmlElement(self):
        element=ElementTree.Element('exclude')
        element.set('date',self.date.strftime('%Y-%m-%dT%H:%M:%S'))
        if self.comment:
            element.set('comment',self.comment)
        return element

    @staticmethod
    def fromXmlElement( element ):
        assert element.tag == 'exclude'
        date = datetime.strptime(element.get('date'),'%Y-%m-%dT%H:%M:%S')
        comment = element.get('comment','')
        return exclude_obs( date, comment )

    
class model( object ):

    BasicComponents=[offset,velocity,annual,semiannual]

    def __init__(self, station=None, xyz=None, filename=None, loadfile=True ):
        self.refdate = refdate
        self.station = station
        self.xyz=xyz
        self.filename=None
        self.enu_axes=None
        self.dates=None
        self.obs=None
        self.useobs=None
        self.excluded=[]
        self.components=[c(self) for c in self.BasicComponents]
        self.covariance=None
        if filename and loadfile:
            self.load(filename)
        self.filename=filename
        self.saved=self.toXmlString()

    def setStation( self, station, xyz ):
        self.station = station
        self.xyz=xyz
        grs80 = ellipsoid.grs80
        lon,lat,h=grs80.geodetic(xyz)
        self.enu_axes=grs80.enu_axes(lon,lat)

    def toXmlString( self ):
        assert re.match(r'^\w{4}$',self.station)
        code=self.station.upper()
        root=ElementTree.Element('station')
        root.set('code',code)
        root.set('site',code)
        root.set('priority','1')
        root.set('start_date','2000-01-01T00:00:00')
        root.set('end_date','')
        
        spm=ElementTree.Element('coordinate_prediction_model')
        spm.set('ref_date',self.refdate.strftime('%Y-%m-%dT%H:%M:%S'))
        if self.xyz:
            spm.set('X0',str(self.xyz[0]))
            spm.set('Y0',str(self.xyz[1]))
            spm.set('Z0',str(self.xyz[2]))
        components=ElementTree.Element('components')
        for c in self.components:
            components.append(c.xmlElement())
        spm.append(components)
        if self.covariance is not None:
            c=self.covariance
            size=c.shape[0]
            covar=ElementTree.Element('covariance')
            covar.set('size',str(size))
            for i in range(size):
                for j in range(i+1):
                    el=ElementTree.Element('element')
                    el.set('row',str(i))
                    el.set('col',str(j))
                    el.set('value',str(c[i,j]))
                    covar.append(el)
            spm.append(covar)
        if self.excluded:
            excluded=ElementTree.Element('excluded')
            for e in self.excluded:
                excluded.append(e.xmlElement())
            spm.append(excluded)
        root.append(spm)
        root.append(ElementTree.Element('outages'))
        return ElementTree.tostring(root)

    def changed( self ):
        xmlstr=self.toXmlString()
        return xmlstr != self.saved

    def save( self, filename=None ):
        if not self.station or self.xyz is None:
             raise RuntimeError('Cannot save stn_pred_model if station code or coordinates not defined')
        if not filename:
            filename=self.filename
        if not filename:
            raise ValueError('No file name specified for saving station prediction  model')
        if '{code}' in filename:
            filename=filename.replace('{code}',self.station)

        with open(filename,'w') as f:
            xmlstr=self.toXmlString()
            pxmlstr=minidom.parseString(xmlstr).toprettyxml(indent='  ')
            f.write(pxmlstr)
            self.saved = xmlstr

    def load( self, filename ):
        if self.station and '{code}' in filename:
            filename = filename.replace('{code}',self.station)
        if not os.path.exists(filename):
            raise RuntimeError('Station prediction model file '+filename+' does not exist')
        with open(filename) as mf:
            xmlstr=mf.read()
            self.loadFromXml(xmlstr,source="Model file "+filename)

        self.filename=filename
        self.saved=self.toXmlString()

    def loadFromXml( self, xmlstr, source="XML string" ):
        #tree=ElementTree.ElementTree(ElementTree.fromstring(xmlstr))
        #root=tree.getroot()
        root=ElementTree.fromstring(xmlstr)
        if root.tag != 'station_prediction_model':
            raise ValueError(source+' is not a station prediction model file')
        station=root.find('station')
        code = station.get('code','')
        if not code:
            raise ValueError(source+' does not specify a station code')
        if self.station and code != self.station:
            raise ValueError(source+' is not for station '+self.station)
        xyz=[0,0,0]
        for i,axis in enumerate(('X','Y','Z')):
            try:
                xyz[i]=float(station.get(axis,''))
            except:
                raise ValueError(source+' does not define an '+axis+' value')
        self.setStation(code,xyz)

        components=[]
        comproot=root.find('components')
        for c in comproot:
            components.append(base_function.fromXmlElement(self,c))
        self.components=components
        self.addBasicComponents()

        self.covariance=None
        covar=root.find('covariance')
        if covar is not None:
            size=int(covar.get('size'))
            self.covariance=np.zeros((size,size))
            for el in covar:
                if el.tag == 'element':
                    r=int(el.get('row'))
                    c=int(el.get('col'))
                    v=float(el.get('value'))
                    self.covariance[r,c]=v
                    self.covariance[c,r]=v

        self.excluded=[]
        excluded=root.find('excluded')
        if excluded is not None:
            for e in excluded:
                self.excluded.append(exclude_obs.fromXmlElement(e))
        self.setExcludedObs()

    def sortComponents( self ):
        bc=self.BasicComponents
        def keyf(comp):
            ct=type(comp)
            return (
                bc.index(ct) if ct in bc else len(bc),
                int(comp._param[-1]),
                ct.__name__,
                comp._param[-1]
                )
        self.components.sort(key=keyf)

    def addBasicComponents( self ):
        ctypes = [type(c) for c in self.components]
        for btype in self.BasicComponents:
            if btype not in ctypes:
                component=btype(self)
                if btype in (annual,semiannual):
                    component.setEnabled(False)
                self.components.append(component)
        self.sortComponents()

    def addComponent( self, component ):
        if component in self.components:
            return
        if component.model != self:
            raise RuntimeError('Cannot add component for a different model')
        self.components.append(component)
        self.sortComponents()

    def removeComponent( self, component ):
        if component in self.components:
            self.components.remove(component)
            component.model=None


    def setExcludedObs( self ):
        '''
        Used to initially link the excluded obs list to the time series data
        '''
        if self.dates is not None:
            self.useobs=np.array([True]*len(self.dates))
            days=[asday(d) for d in self.dates]
            for e in self.excluded:
                e.index=-1
                for i,d in enumerate(days):
                    if abs(d-e.day) < 0.5:
                        e.index=i
                        self.useobs[i]=False
                        break

if __name__ == '__main__':
    import argparse

    parser=argparse.ArgumentParser('Convert station prediction model XML to station XML')
    parser.add_argument('codes',nargs='+',help='Codes of stations to analyse')
    parser.add_argument('-m','--model-dir',default='models',help='Base directory for models')
    parser.add_argument('-s','--station-dir',default='stations',help='Base directory for stations')

    args=parser.parse_args()
    codes = [c.upper() for c in args.codes]

    model_file=args.model_dir+'/{code}_spm.xml'
    station_file=args.station_dir+'/{code}.xml'

    assert os.path.isdir(args.model_dir), args.model_dir + ' is not a directory'
    assert os.path.isdir(args.station_dir), args.station_dir + ' is not a directory'

    if codes == ['*']:
        codes=[]
        for f in os.listdir(args.model_dir):
            match = re.match(r'^(\w{4})_spm\.xml',f)
            if match:
                codes.append(match.group(1))

    for code in codes:
        print "Converting",code
        m=model(station=code,filename=model_file)
        m.save(filename=station_file)


