# positionz_coordinate_models

This repository contains the station coordinate models representing the trajectory of the PositioNZ stations.

These models are updated periodically (6 months to 1 year) by manually calculating using the python spm_editor
software at https://github.com/linz/python-linz-stationcoordmodel.  The models are based on the ITRF2008 coordinate
time series derived from the LINZ daily processing of the PositioNZ station RINEX data.

The models are used to calculate the coordinates of the stations in the LINZ online GNSS post-processing service
at http://www.linz.govt.nz/positionzpp.
