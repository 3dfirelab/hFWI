from datetime import datetime, timedelta
from dataHandler import load_patch_for_date
from time import time
from os import path, makedirs
import numpy as np
import sys
import pdb 

from fwi import Indices

latitude_arome = 46.45


PATCH_SIZE=23
PATCH_IDX=int(sys.argv[1])

start    = datetime.strptime('2023-12-31T00:00:00Z', "%Y-%m-%dT%H:%M:%SZ")
endspin  = datetime.strptime('2024-03-31T00:00:00Z', "%Y-%m-%dT%H:%M:%SZ")
stop     = datetime.strptime('2024-04-04T00:00:00Z', "%Y-%m-%dT%H:%M:%SZ")
date = start

wind, rh, temp, rain = load_patch_for_date(date, PATCH_IDX)
indices = Indices(wind.shape, latitude_arome, PATCH_IDX)

while date <= stop:
    print(date.strftime('%Y-%m-%dT%H:%M:%SZ'))
    try:
        wind, rh, temp, rain = load_patch_for_date(date, PATCH_IDX)
        #print("       ", end=' ')
    except:
        pass
        #print("missing", end=' ')
        
    time2dateStop = (stop-date).total_seconds()
    IntMode = 'hourly'
    if time2dateStop > 48 * 3600:
        IntMode = 'daily'

    indices.integrate(date, IntMode, wind, rh, temp, rain)

    #print(date.hour, IntMode,'{:.2f}'.format(indices.fwi.mean()))
        
    if date.hour == 11 or date >= endspin:
        out_dir = f"/home/fire/output-rp/{date.strftime('%Y-%m-%dT%H:%M:%SZ')}"
        makedirs(out_dir, exist_ok=True)
            
        out_file = path.join(out_dir, f"patch-{PATCH_IDX}.npz")
        np.savez(out_file,
                ffmc = indices.ffmc,
                dmc = indices.dmc,
                dc = indices.dc,
                isi = indices.isi,
                bui = indices.bui,
                fwi = indices.fwi
                )

    
    date = date + timedelta(hours=1)



