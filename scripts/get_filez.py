import numpy as np
import subprocess
import astropy as ap
from astropy.io import fits
from astropy import table as t
from astropy.table import Table
import os
import time

os.chdir("/Volumes/Nitya/data/mpl8")

location = '/Volumes/Nitya/data/mpl8'

def bash_command(url):
    p = subprocess.Popen(['wget','--user', 'sdss', '--password',
                                '2.5-meters','-P',location,  url],
            stdout=subprocess.PIPE)
    out, err = p.communicate()
    if err:
        print(err)

drpall = t.Table.read('drpall-v2_5_3.fits')
drpall.add_index('plateifu')
index1 = np.where(drpall['srvymode']=='MaNGA dither')[0]
drpall1 = drpall[index1]
index2 = np.where(drpall1['z']>0)[0]
drpall = drpall1[index2]

plateifu = drpall['plateifu']

earl = ['https://data.sdss.org/sas/mangawork/manga/spectro/redux/MPL-8/' +
        str(drpall['plate'][i]) + '/stack/manga-' +  plateifu[i]  +
        '-LOGCUBE.fits.gz' for i in range(4660, 6491)]

#subprocess.Popen(['wget', '--user', 'sdss', '--password',
#'2.5-meters', '-P', location,  'https://data.sdss.org/sas/mangawork/manga/spectro/redux/MPL-7/
#8935/stack/manga-8935-12701-LOGCUBE.fits.gz'], stdout = subprocess.PIPE)

#earl = np.asarray(earl)
for i in range(len(earl)):
        bash_command(earl[i])
        if i!=0 and i%100==0:
            print(i , earl[i], 'galaxy donezos but we be waiting' )
            time.sleep(120)
            print('ayyeee wait over')
