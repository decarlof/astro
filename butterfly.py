import json, urllib, requests, numpy as np, matplotlib.pylab as plt, matplotlib.ticker as mtick
from datetime import datetime as dt_obj
import matplotlib.dates as mdates
import matplotlib.colors as mcol
from matplotlib.dates import  * 
from astropy.io import fits
import math
import wget
import os 


def parse_tai_string(tstr, datetime=True):
    year   = int(tstr[:4])
    month  = int(tstr[5:7])
    day    = int(tstr[8:10])
    hour   = int(tstr[11:13])
    minute = int(tstr[14:16])
    if datetime: return dt_obj(year, month, day, hour, minute)
    else: return year, month, day, hour, minute

def download_info(url, query=''):

    # get the file names for the data files from the beginning of the mission to the present moment
    response = requests.get(url + query)
    data     = response.json()
    t_last   = data['keywords'][0]['values'][0]
    
    t_rec = None
    mission_data = None
    local_filename = None
    if query != '':
        t_rec    = data['keywords'][0]['values']    # create a list of the keyword t_rec
        remote_filename = data['segments'][0]['values']    # create a list of the remote_filename

        # print some statistics
        n_elements = len(remote_filename)
        print("There are "       + str(n_elements) + " files.")
        print("first record is " + t_rec[0]        + " file=" + remote_filename[0])
        print("last record is "  + t_rec[-1]       + " file=" + remote_filename[-1])

        local_filename=[]
        for i in range(len(remote_filename)):
            local_filename.append(remote_filename[i].replace('/','_'))

        mission_data = {}
        for key in range(len(remote_filename)):
            mission_data[key] = [t_rec[key], remote_filename[key], local_filename[key]]

    return t_last, mission_data

def dowload_the_polar_field_data(mission_data, remote=False):

    t_rec = []
    remote_filename = []
    local_filename  = []
    for i in range(len(mission_data)):
       t_rec.append(mission_data[i][0])
       remote_filename.append(mission_data[i][1])
       local_filename.append(mission_data[i][2])

    if remote == True:
        for i in range(len(mission_data)):
            url = "http://jsoc.stanford.edu" + remote_filename[i]
            fname = wget.download(url)
            os.rename(fname, local_filename[i])
            print('\n'+str(remote_filename[i]))

def gather_the_polar_field_data(mission_data, local=True):

    t_rec = []
    remote_filename = []
    local_filename  = []
    for i in range(len(mission_data)):
        t_rec.append(mission_data[i][0])
        remote_filename.append(mission_data[i][1])
        local_filename.append(mission_data[i][2])

    # format t_rec into a numpy array
    t_rec = np.array(t_rec, dtype='S16') # this dtype is ok to set as the format of T_REC will never change

    # convert t_rec into a datetime object
    x = np.array([parse_tai_string(t_rec[i], datetime=True) for i in range(t_rec.size)])

    # convert the datetime object into a floating point number
    xmdates = mdates.date2num(x)
    timespan = xmdates[-1] - xmdates[0]
    nsteps = int(timespan) + 1    
    
    # use the astropy library to open the FITS data
    mf_br = np.ndarray([180, nsteps])
    mf_br.fill(np.nan)
    for i in range(len(remote_filename)):
        itime = int(xmdates[i] - xmdates[0])
        if local:
            hdu = fits.open(local_filename[i])
        else:
            url = "http://jsoc.stanford.edu"+remote_filename[i]
            hdu = fits.open(url)
        mf_br[:, itime] = hdu[0].data
        
    # find the pixels with values greater than 200 G and set them equal to NaN
    qq = np.where(abs(mf_br) > 200) # this throws an error because there are nans in the array but it still works
    print("There are", len(qq[0]), "pixels greater than 200 G.")
    for i in range(len(qq[0])):
        mf_br[qq[0][i], qq[1][i]] = np.nan
    qq = np.where(abs(mf_br) > 200) # this throws an error because there are nans in the array but it still works
    print("Now there are", len(qq[0]), "pixels greater than 200 G.")

    return mf_br, xmdates, nsteps

def differential_rotation_period(latitude):
    k = 3.6e3 * 24. * 1e-6 * (180./math.pi)
    a = 2.7139 * k
    b = -0.405 * k
    c = -0.422 * k
    theta = (90. - latitude) * (math.pi/180.)
    cos_theta = math.cos(theta) 
    omega = a  +  b * cos_theta * cos_theta  +  c * cos_theta * cos_theta * cos_theta * cos_theta
    period = 360./omega
    return period

def process_the_polar_field_data(mf_br, xmdates, nsteps, n_elements):
    mf_br_averaged = np.ndarray([180, nsteps])
    mf_br_averaged.fill(np.nan)
    for j in range(n_elements):
        jtime = int(xmdates[j] - xmdates[0])
        for i in range(180):
            latitude = i-90.
            period = differential_rotation_period(latitude)
            period_even = math.ceil(period/2.) * 2 # round the period to the nearest even number
            half_period_even = int(math.ceil(period/2.) * 2)/2
            mf_br_averaged[i, jtime] = np.nanmean(mf_br[i, int(jtime-half_period_even):int(jtime + half_period_even)])
        if j % 365 == 0:
            print("j=", str(jtime), " ", np.nanmean(mf_br_averaged[:, jtime]))    
    return mf_br_averaged, xmdates, nsteps

def plot_butterfly(mf_br_averaged, xmdates, nsteps):
    
    # create the figure
    fig, ax = plt.subplots(1, 1)

    # define some style elements
    marker_style = dict(linestyle='-', linewidth=4, fillstyle='full')
    text_style = {'color':'black', 'weight':'normal', 'size': 16}

    # create a customized color table
    colors = np.loadtxt('cstretch.txt')
    colormap_customized = mcol.ListedColormap(colors)

    # plot the array
    im = plt.imshow(mf_br_averaged, cmap=colormap_customized, vmax=6, vmin=-6, 
                    extent=[xmdates[0], xmdates[-1], -90, 90], origin='lower')

    # format the x-axis with universal time
    ax.xaxis_date()
    locator = AutoDateLocator()
    locator.intervald[MONTHLY] = [6]
    formatter = DateFormatter('%Y')
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    # label the axes and the plot
    ax.set_xlabel('time', labelpad=10, fontdict = text_style)
    ax.set_ylabel('latitude (deg)', fontdict=text_style)
    ax.set_title('zonally averaged radial component of magnetic field', fontdict = text_style)
    ax.set_aspect(5)
    fig.set_figwidth(nsteps/365.25)
    
    return fig

def main():

    # get the most recent timestamp
    url      = "http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info?ds=hmi.meanpf_720s[$]&op=rs_list&key=T_REC"
    t_last, _ = download_info(url)

    # get the file names for the data files from the beginning of the mission to the present moment (t_last)
    url      = "http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info?"
    query    = "ds=hmi.meanpf_720s[2010.05.01_00_TAI-" + t_last + "@1d]&op=rs_list&key=T_REC&seg=mf_br"

    _, new_mission_data = download_info(url, query)

    dowload_the_polar_field_data(new_mission_data, remote=True)

    hmidata, hmidates, hmisteps = gather_the_polar_field_data(new_mission_data)
    hmidata, hmidates, hmisteps = process_the_polar_field_data(hmidata, hmidates, hmisteps, len(new_mission_data))
    fig      = plot_butterfly(hmidata, hmidates, hmisteps)

    # get the file names for the data files from the SoHO/MDI data at a sampling of once per day
    url = "http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info?"
    query = "ds=mdi.meanpf_96m[1996.05.01_00_TAI-2010.05.01_00_TAI@1d]&op=rs_list&key=T_REC&seg=mf_br"

    _, old_mission_data = download_info(url, query)
    
    dowload_the_polar_field_data(old_mission_data, remote=True)

    mdidata, mdidates, mdisteps = gather_the_polar_field_data(old_mission_data)
    mdidata, mdidates, mdisteps = process_the_polar_field_data(mdidata, mdidates, mdisteps, len(old_mission_data))

    fig = plot_butterfly(mdidata, mdidates, mdisteps)

    hmidata_stronger = 1.2  *  hmidata
    totaldata = np.ndarray([180, mdidata.shape[1] + hmidata_stronger.shape[1]])
    totaldata[:, 0:mdidata.shape[1]] = mdidata
    totaldata[:, mdidata.shape[1]:mdidata.shape[1] + hmidata_stronger.shape[1]] = hmidata_stronger
    print("The total data array has the following dimensions:", totaldata.shape, ".")

    totalsteps = hmisteps  +  mdisteps
    print("There are", totalsteps, "total time steps.")

    totaldates = np.arange(mdidates[0], hmidates[-1])

    fig = plot_butterfly(totaldata, totaldates, totalsteps)

    fig.savefig('butterfly.png', dpi=300)

if __name__ == "__main__":
    main()