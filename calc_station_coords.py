
defdir='/home/ccrook/projects/deformation/models/published'
import sys
sys.path.append(defdir+'/tools')

import os
import re
import argparse
from datetime import datetime
from stn_pred_model import model as spm
from LINZ.DeformationModel import Model as DefModel
from LINZ.Geodetic.ellipsoid import grs80
from LINZ.Geodetic.ITRF_transformation import transformation


parser=argparse.ArgumentParser('Script to compare NZGD2000 coordinates with station time series')
parser.add_argument('gdb_file',default='gdb_coords.txt',help="Name of gdb file input file")
parser.add_argument('check_file',nargs='?',default='coord_differences.csv',help="Output csv file of coordinate comparisons")
parser.add_argument('update_csv_file',nargs='?',help="Name of gdb coordinate upload file (default none)")

args=parser.parse_args()

gdbfile=args.gdb_file
chkfile=args.check_file
updfile=args.update_csv_file

calcDate=datetime(2014,1,1)
itrf_tfm=transformation(from_itrf='ITRF2008',to_itrf='ITRF96')
defmodel=DefModel.Model(defdir+'/model')

gdbcrds={}
with open(gdbfile,'r') as gdbf:
    l=gdbf.readline()
    l=l.lower().replace("\"","").replace(","," ")
    if l.split() == 'code gdb_lon gdb_lat gdb_hgt'.split():
        crdorder=(1,2,3)
    elif l.split() == 'geodeticcode lat lon ellhgt'.split():
        crdorder=(2,1,3)
    else:
        raise RuntimeError("Invalid fields in "+gdbfile)
    for l in gdbf:
        l=l.lower().replace("\"","").replace(","," ")
        try:
            parts=l.split()
            crds=[float(parts[x]) for x in crdorder]
            gdbcrds[parts[0].upper()]=crds
        except:
            pass

csvu=None
if updfile:
    csvu=open(updfile,'w')

with open(chkfile,'w') as csv:
    csv.write(','.join((
        'code',
        'itrf2008_X','itrf2008_Y','itrf2008_Z',
        'itrf2008_lon','itrf2008_lat','itrf2008_hgt',
        'itrf96_lon','itrf96_lat','itrf96_hgt',
        'nzgd2000_lon','nzgd2000_lat','nzgd2000_hgt',
        'gdb_lon','gdb_lat','gdb_hgt',
        'e_diff','n_diff','h_diff',
    )))
    csv.write('\n');
    if csvu:
        csvu.write("CODE,LAT,LON,EHGT,ROBG,COSY,DATE\n");


    codes=[]
    for f in os.listdir('stations'):
        m=re.match(r'^(\w{4}).xml$',f)
        if not m:
            continue
        codes.append(m.group(1))

    for code in sorted(codes):
        f=code+'.xml'
        if code not in gdbcrds:
            continue

        print "Processing",code

        try:
            m=spm(filename='stations/'+f)
            # Don't want annual and semi-annual components...
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

            ucode=code.upper()
            if ucode in gdbcrds:
                gcrds=gdbcrds[ucode] 
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*gcrds))
                dedln,dndlt=grs80.metres_per_degree(*gcrds)
                edif=(gcrds[0]-llhnz2k[0])*dedln
                ndif=(gcrds[1]-llhnz2k[1])*dndlt
                hdif=gcrds[2]-llhnz2k[2]
                csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(edif,ndif,hdif))
            else:
                csv.write(',,,,,')
            csv.write("\n")
            if csvu and ucode in gdbcrds:
                csvu.write('"{0}",{1:.9f},{2:.9f},{3:.4f},"B10","NZGD2000","2000.01.01"\n'.format(
                    code.upper(),*llhnz2k))
        except:
            print sys.exc_info()[1]
