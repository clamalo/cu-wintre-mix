import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import pyart
from mpl_toolkits.axes_grid1 import make_axes_locatable


directories = ['CAND/20220222_IOP05/NetCDF/VER/low/','USD/20220222_IOP05/NetCDF/VER/high/','USD/20220222_IOP05/NetCDF/VER/low/']
# # '/Volumes/Data/2022_WINTRE-MIX/radar/COW/20220222_IOP05/NetCDF/VER/high/20220222/'
# # '/Volumes/Data/2022_WINTRE-MIX/radar/COW/20220222_IOP05/NetCDF/VER/high/20220223/'
# # '/Volumes/Data/2022_WINTRE-MIX/radar/COW/20220222_IOP05/NetCDF/VER/low/20220222/'
# # '/Volumes/Data/2022_WINTRE-MIX/radar/COW/20220222_IOP05/NetCDF/VER/low/20220223/'
# '/Volumes/Data/2022_WINTRE-MIX/radar/USD/20220222_IOP05/NetCDF/VER/high/'
# '/Volumes/Data/2022_WINTRE-MIX/radar/USD/20220222_IOP05/NetCDF/VER/low/'


# directory = '/Volumes/Data/2022_WINTRE-MIX/radar/CAND/20220222_IOP05/NetCDF/VER/low/'

for directory in directories:
    file_list = os.listdir('/Volumes/Data/2022_WINTRE-MIX/radar/'+directory)
    print(len(file_list))

    parameters = ['DBZHC_F', 'ZDRM', 'KDP', 'RHOHV','SNRHC', 'WIDTH', 'VEL']
    for parameter in parameters:

        master_values = []
        n = 0
        for file in file_list:
            # file = file_list[n]

            n+=1
            print(n)

            print(file)

            if file == 'cfrad.20220222_214810.628_to_20220222_214822.699_DOW7low_VER.nc' or file == 'low.zip' or '.nc' not in file:
                pass
            else:
                ds = xr.load_dataset('/Volumes/Data/2022_WINTRE-MIX/radar/'+directory+file,engine='netcdf4')



                values = []
                y_values = []

                for r in range(len(ds.range.values)):
                    range_ds = ds.isel(range=r)

                    # print(len(range_ds.DBZHC_F.values))
                    time_mean = np.mean(range_ds[parameter].values)
                    values.append(time_mean)
                    y_values.append(r*74.94810968637466)

                master_values.append(values)


        #create numpy array from 0 to 77
        times = np.arange(0,len(master_values),1)
        print(times)
        master_ds = xr.Dataset({
            parameter : xr.DataArray(
                        data   = master_values,   # enter data here
                        dims   = ['time','height'],
                        coords = {'time': times},
                        attrs  = {
                            '_FillValue': -999.9,
                            'units'     : 'N/A'
                            }
                        ),
                    },
                attrs = {'example_attr': 'this is a global attribute'}
            )

        master_ds.to_netcdf('/Users/clamalo/documents/'+parameter+'master_ds.nc')

        ds = xr.load_dataset('/Users/clamalo/documents/'+parameter+'master_ds.nc')
        ds['height'] = ds['height']*74.94810968637466/1000
        ds = ds.isel(height=slice(0,77))
        ds = ds.transpose('height','time')
        print(ds)

        #set figsize to be a square
        # plt.figure(figsize=(10,10))
        cf = plt.pcolormesh(ds.time.values, ds.height.values, ds[parameter].values)
        plt.title('IOP05 - '+parameter+' (89 deg scan)')
        plt.xlabel('Time')
        plt.ylabel('Height (m)')
        plt.colorbar(cf)
        plt.savefig('/Users/clamalo/documents/'+parameter+'.png', dpi=1000)
        plt.clf()



    parameters = ['DBZHC_F', 'ZDRM', 'KDP', 'RHOHV','SNRHC', 'WIDTH', 'VEL']
    master_ds = xr.load_dataset('/Users/clamalo/documents/DBZHC_Fmaster_ds.nc')
    for parameter in parameters:
        ds = xr.load_dataset('/Users/clamalo/documents/'+parameter+'master_ds.nc')
        master_ds[parameter] = ds[parameter]

    new_directory = directory.replace('/','_')
    master_ds.to_netcdf('/Users/clamalo/documents/'+new_directory+'IOP05.nc')

    master_ds['height'] = master_ds['height']*74.94810968637466
    master_ds = master_ds.isel(height=slice(0,77))
    master_ds = master_ds.transpose('height','time')

    parameters = ['DBZHC_F', 'ZDRM', 'KDP', 'RHOHV','SNRHC', 'WIDTH', 'VEL']

    # master_ds['time'] = master_ds['time']*6

    valid_times = []
    m = 0
    for file in file_list:
        if '.nc' in file:
            valid_time = file.split('.')[1]
            year = valid_time[0:4]
            month = valid_time[4:6]
            day = valid_time[6:8]
            hour = valid_time[9:11]
            minute = valid_time[11:13]
            # valid_time = year+'-'+month+'-'+day+' '+hour+':'+minute
            valid_time = hour+':'+minute
            if m%5 == 0:
                valid_times.append(valid_time)
            m+=6

    times = master_ds.time.values
    heights = master_ds.height.values
    DBZHC = master_ds.DBZHC_F.values
    ZDRM = master_ds.ZDRM.values
    KDP = master_ds.KDP.values
    RHOHV = master_ds.RHOHV.values
    SNRHC = master_ds.SNRHC.values
    WIDTH = master_ds.WIDTH.values
    VEL = master_ds.VEL.values
    # cf = plt.pcolormesh(ds.time.values, ds.height.values, ds[parameter].values)

    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    axs = [ax1,ax2,ax3]

    for ax in axs:
        if ax == ax1:
            im = ax.pcolormesh(times, heights, SNRHC, vmin=-20, vmax=30, cmap= plt.get_cmap('pyart_NWSRef'))
        elif ax == ax2:
            im = ax.pcolormesh(times, heights, WIDTH, vmin=-1, vmax=2, cmap = plt.get_cmap('pyart_NWS_SPW'))
        elif ax == ax3:
            im = ax.pcolormesh(times, heights, VEL, vmin=-15, vmax=15, cmap = plt.get_cmap('pyart_NWSVel'))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        if ax != ax3:
            ax.set_xticks(np.arange(0, len(times)+1, 5.0))
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xticks(np.arange(0, len(times)+1, 5.0),labels=valid_times, rotation = 45, fontsize=12)
            # ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        ax.set_xlim(0,times[-1])
        ax.set_ylim(0,heights[-1])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        if ax == ax1:
            cbar.set_label('Signal-to-noise ratio (dB)', rotation=90, labelpad=8, fontsize=15)
            cbar.set_ticks([-20, -10, 0, 10, 20, 30])
        elif ax == ax2:
            cbar.set_label('Spectrum Width (m/s)', rotation=90, labelpad=8,fontsize=15)
            cbar.set_ticks([-1, -0.5, 0, 0.5, 1, 1.5, 2])
        elif ax == ax3:
            cbar.set_label('Velocity (m/s)', rotation=90, labelpad=8,fontsize=15)
            cbar.set_ticks([-15,-10,-5,0,5,10,15])
        cbar.ax.tick_params(labelsize=15)
        ax.grid(color='black',linestyle='-.')
    fig.text(0.5, 0.04, 'Time', ha='center', fontsize=15)
    fig.text(0.05, 0.5, 'Height AGL (m)', va='center', rotation='vertical', fontsize=15)
    new_directory = directory.replace('/','_')
    plt.savefig('/Users/clamalo/documents/'+new_directory+'_3.png', dpi=300, bbox_inches='tight')
    plt.clf()

    fig = plt.figure(figsize=(12, 16))
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)

    axs = [ax1,ax2,ax3,ax4]

    for ax in axs:
        if ax == ax1:
            im = ax.pcolormesh(times, heights, DBZHC, vmin=-20, vmax=30, cmap = 'pyart_NWSRef')
        elif ax == ax2:
            im = ax.pcolormesh(times, heights, ZDRM, vmin=-1, vmax=2, cmap = 'pyart_NWSRef')
        elif ax == ax3:
            im = ax.pcolormesh(times, heights, KDP, vmin=0, vmax=2, cmap='plasma')
        elif ax == ax4:
            im = ax.pcolormesh(times, heights, RHOHV, vmin=0.8, vmax=1, cmap = plt.get_cmap('pyart_RdYlBu11b_r'))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        if ax != ax4:
            ax.set_xticks(np.arange(0, len(times)+1, 5.0))
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xticks(np.arange(0, len(times)+1, 5.0),labels=valid_times, rotation = 45, fontsize=12)
            # ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        ax.set_xlim(0,times[-1])
        ax.set_ylim(0,heights[-1])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        if ax == ax1:
            cbar.set_label('Z (dBZ)', rotation=90, labelpad=8, fontsize=15)
            cbar.set_ticks([-20, -10, 0, 10, 20, 30])
        elif ax == ax2:
            cbar.set_label('ZDR (dB)', rotation=90, labelpad=8,fontsize=15)
            cbar.set_ticks([-1, -0.5, 0, 0.5, 1, 1.5, 2])
        elif ax == ax3:
            cbar.set_label('KDP (deg/km)', rotation=90, labelpad=8,fontsize=15)
            cbar.set_ticks([0, 0.3, 0.8, 1.2, 1.6, 2])
        elif ax == ax4:
            cbar.set_label('RHOHV', rotation=90, labelpad=8,fontsize=15)
            cbar.set_ticks([0.8, 0.85, 0.9, 0.95, 1])
        cbar.ax.tick_params(labelsize=15)
        ax.grid(color='black',linestyle='-.')
    fig.text(0.5, 0.06, 'Time', ha='center', fontsize=15)
    fig.text(0.05, 0.5, 'Height AGL (m)', va='center', rotation='vertical', fontsize=15)
    new_directory = directory.replace('/','_')
    plt.savefig('/Users/clamalo/documents/'+new_directory+'_4.png', dpi=300, bbox_inches='tight')
    plt.clf()