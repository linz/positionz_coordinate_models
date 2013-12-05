import stn_pred_model
import os
import os.path
import re
import sys

model_dir='stations'
timeseries_file='timeseries/{code}_igs08_xyz.dat'


for f in sorted(os.listdir(model_dir)):
    m=re.match('^(\w{4})\.xml$',f)
    if not m:
        print "Skipping file",f
        continue
    code=m.group(1)
    tsf = timeseries_file.replace('{code}',code)
    if not os.path.exists(tsf):
        print "No time series data for",code
        continue
    print "Updating code",code
    try:
        model=stn_pred_model.model(filename=os.path.join(model_dir,f))
        model.loadTimeSeries(tsf)
        model.save(updateAvailability=True)
    except:
        print "Error:",str(sys.exc_info()[1])
