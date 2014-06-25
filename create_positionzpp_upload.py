#!/usr/bin/python

import argparse
import re
import os
import os.path
import md5
import zipfile

parser=argparse.ArgumentParser(description="Create a PositioNZ-PP station prediction model upload file")
parser.add_argument('zip_file',help="Name of zip file")
parser.add_argument('spm_xml_file',nargs='*',help='Station prediction model XML files to upload CCCC.xml')
parser.add_argument('-r','--remove-file',nargs='*', help='The station code of files to remove from positionzpp')

args=parser.parse_args()
upload_key='GNSSSPM'

codes=[]
files=[]
spmhash=''
for sfile in args.spm_xml_file:
    if not os.path.exists(sfile):
        print sfile,'does not exist!'
        continue
    filename=os.path.basename(sfile)
    if not re.match(r'^\w{4}\.xml$',filename):
        print sfile,' is not a valid filename for a station prediction model'
        continue
    code=filename[0:4].upper()
    codes.append(code)
    with open(sfile) as sf:
        data=sf.read()
        m=md5.new()
        m.update(upload_key)
        m.update(data)
        spmhash=spmhash+code+".xml "+m.hexdigest()+"\n"
    files.append({'file':code+'.xml','data':data})

remove_files=args.remove_file or []
for rcode in remove_files:
    if not re.match(r'^\w{4}$',rcode):
        print "Invalid code".rcode," for removal"
        continue
    rcode=rcode.upper()
    if rcode in codes:
        print "Cannot remove {0} as already used".format(rcode)
        continue
    m=md5.new()
    m.update(upload_key)
    m.update("REMOVE")
    m.update(rcode)
    m.update('.xml')
    spmhash=spmhash+rcode+".xml "+m.hexdigest()+"\n"

zip=zipfile.ZipFile(args.zip_file,'w')
zip.writestr('spm.hash',spmhash)
for f in files:
    zip.writestr(f['file'],f['data'])
zip.close()







