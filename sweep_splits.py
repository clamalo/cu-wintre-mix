import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import pyart

directory = '/Volumes/Data/2022_WINTRE-MIX/radar/CAND/20220222_IOP05/NetCDF/VOL/SURtest_0.5deg/low/'

file_list = os.listdir(directory)
print(len(file_list))

current_cmap = plt.get_cmap('pyart_NWSRef')
current_min = -20
current_max = 30

tilts = [0.41,0.8,1.2,1.61,1.99,2.8,3.59]
for t in range(len(file_list)):
    os.mkdir('/Users/clamalo/documents/radar/'+str(t))
    # for tilt in tilts:
    #     os.mkdir('/Users/clamalo/documents/radar/'+str(t)+'/'+str(tilt))


t = 8
for file in file_list:

    file = file_list[t]
    if 'COW' in directory+file:
        max_range = 1181
    else:
        max_range = 1000

    print(t)

    ds = xr.load_dataset(directory+file)

    starts = list(ds.sweep_start_ray_index.values[:-1])
    ends = list(ds.sweep_end_ray_index.values[:-1])

    starts.reverse()
    ends.reverse()

    #slice of time
    ds = ds.isel(time=slice(0,int(ends[0])+1))

    tilts = [0.41,0.8,1.2,1.61,1.99,2.8,3.59]

    for n in range(len(tilts)):
        print((int(starts[n])),' ',(int(ends[n])+1))
        new_ds = ds.isel(time=slice((int(starts[n])),(int(ends[n])+1)))
        # new_ds.to_netcdf('/Users/clamalo/downloads/'+str(tilts[n])+'_radar.nc')
        values = new_ds.DBZHC.values
        
        # Using linspace so that the endpoint of 360 is included...
        print(ends[n]+1-starts[n])
        if ends[n]+1-starts[n] == 719:
            num = 719
        else:
            num = 720
        azimuths = np.radians(np.linspace(0, 360, num))

        zeniths = np.arange(0, max_range, 1)

        r, theta = np.meshgrid(zeniths, azimuths)

        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        levels = np.linspace(-20, 30, 50)
        cf = ax.contourf(theta, r, values, cmap = current_cmap, vmin = current_min, vmax = current_max, levels = levels)
        cbar = plt.colorbar(cf)

        # plt.savefig('/Users/clamalo/documents/radar/'+str(t)+'/'+str(tilts[n])+'/radar.png', dpi=300)
        plt.savefig('/Users/clamalo/documents/radar/'+str(t)+'/'+str(tilts[n])+'_radar.png', dpi=300)
        plt.clf()

        if 'COW' in directory+file:
            ranges = [67,134,201,268,335,402,469,536,603,670,737,804,871,938,1005,1072,1139]
        else:
            ranges = [67,134,201,268,335,402,469,536,603,670,737,804,871,938]
        k = 0
        for beam_range in ranges:
            k = k+1
            new_new_ds = new_ds.isel(range=beam_range)
            range_values = []
            if ends[n]+1-starts[n] == 719:
                for p in range(719):
                    range_values.append(new_new_ds.DBZHC_F.values[p])
            else:
                for p in range(720):
                    range_values.append(new_new_ds.DBZHC_F.values[p])
            #line thickness
            plt.plot(range_values, linewidth=0.5)
        plt.xlim(0,720)
        plt.savefig('/Users/clamalo/documents/radar/'+str(t)+'/'+str(tilts[n])+'_all_tilts.png', dpi=300)
        plt.clf()
    t+=1