
defdir='/home/ccrook/projects/deformation/models/published'
import sys
sys.path.append(defdir+'/tools')

import os
import re
from datetime import datetime
from stn_pred_model import model as spm
from LINZ.DeformationModel import Model as DefModel
from LINZ.Geodetic.ellipsoid import grs80
from LINZ.Geodetic.ITRF_transformation import transformation


file='stations/AUCK.xml'
calcDate=datetime(2014,1,1)
itrf_tfm=transformation(from_itrf='ITRF2008',to_itrf='ITRF96')
defmodel=DefModel.Model(defdir+'/model')


with open('positionz_coordinates.csv','w') as csv:
    csv.write(','.join((
        'code',
        'itrf2008_X','itrf2008_Y','itrf2008_Z',
        'itrf2008_lon','itrf2008_lat','itrf2008_hgt',
        'itrf96_lon','itrf96_lat','itrf96_hgt',
        'nzgd2000_lon','nzgd2000_lat','nzgd2000_hgt'
    )))
    csv.write('\n');

    for f in os.listdir('stations'):
        m=re.match(r'^(\w{4}).xml$',f)
        if not m:
            continue
        code=m.group(1)

        print "Processing",code

        try:
            m=spm(filename='stations/'+f)
            for c in m.components:
                if 'annual' in c.componentType():
                    c.setEnabled(False)

            xyz08=m.calc(calcDate,enu=False)
            llh08=grs80.geodetic(xyz08)
            llh96=itrf_tfm.transformLonLat(llh08[0],llh08[1],llh08[2],calcDate)
            if llh96[0] < 0:
                llh96[0] += 360.0
            llhnz2k=defmodel.applyTo(llh96,date=calcDate,subtract=True)
            if llh96[0] > 180:
                llh96[0] -= 360.0
            if llhnz2k[0] > 180:
                llhnz2k[0] -= 360.0

            csv.write('"{0}"'.format(code))
            csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(*xyz08))
            csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh08))
            csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh96))
            csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llhnz2k))
            csv.write("\n")
        except:
            print sys.exc_info()[1]





