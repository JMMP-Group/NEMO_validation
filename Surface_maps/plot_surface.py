from config import config

cfg = config() # initialise variables in python

import coast
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import cmocean

class plot_surface(object):
    def __init__(self):

        # paths
        clim_path = "/surface_maps/surface_state_climatology.nc"
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_dat = cfg.dn_out + clim_path
        en4_nc = "/surface_maps/en4_gridded_surface_climatology.nc"
        self.en4_grid = cfg.dn_out + en4_nc

    def set_cbar(self, ax, p, txt):
        pos = ax.get_position()
        cbar_ax = fig.add_axes([pos.x0, 0.08, 
                                pos.x1 - pos.x0, 0.02])
        cbar = fig.colorbar(p, cax=cbar_ax, orientation='horizontal')
        cbar.ax.text(0.5, -3.0, txt, fontsize=8,
                     rotation=0, transform=cbar.ax.transAxes,
                     va='top', ha='center')

    def render(self, ax, da, crange):
        if da.name == 'temperature':
            cmap = cmocean.cm.thermal
        if da.name == 'salinity':
            cmap = cmocean.cm.haline

        p = ax.pcolor(da.longitude_bins, da.latitude_bins, da,
                      vmin=crange[0], vmax=crange[1], cmap=cmap,
                      shading='nearest', transform=self.plt_proj)
        ax.coastlines()
        ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
        return p

    def plot_surface_state_climatology(self, scalar="temperature"):
        """ plot surface state climatology over seasonal means """

        # get data and alias
        c_ds = xr.open_dataset(self.fn_dat)

        tmin = np.floor(c_ds.temperature.min().values)
        tmin = 5
        tmax = np.ceil(c_ds.temperature.max().values)
        smin = 20
        smax = np.ceil(c_ds.salinity.max().values)

        proj=ccrs.AlbersEqualArea()
        self.plt_proj=ccrs.PlateCarree()
        proj_dict = {"projection": self.plt_proj}
        fig, axs = plt.subplots(4, 2, figsize=(5.5,5.5), subplot_kw=proj_dict)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15,
                            hspace=0.05, wspace=0.05)

        seasons = ['DJF','MAM','JJA','SON']
        for i, season in enumerate(seasons):
            da = c_ds.sel(season=season).load()
            print (da)
            pt = self.render(axs[i,0], da.temperature, crange=(tmin,tmax))
            ps = self.render(axs[i,1], da.salinity, crange=(smin,smax))

        self.set_cbar(axs[-1,0], pt, r'Surface Temperature ($^{\circ}$C)')
        self.set_cbar(axs[-1,1], ps, r'Surface Salinity (-)')

        plt.savefig('FIGS/surface_climatology.png', dpi=600)


    def plot_surface_state_climatology_all(self):
        """ plot all climatology for temperature and salinity """

        self.plot_suface_state_climatology("temperature")
        self.plot_suface_state_climatology("salinity")

    def plot_gridded_surface_en4(self):
        """ plot en4 data on regular lat lon grid """

        # get data
        c_ds = xr.open_dataset(self.en4_grid)
        c_ds = c_ds.rename(dict(potential_temperature="temperature",
                                practical_salinity="salinity"))
        print (c_ds)

        tmin = np.floor(c_ds.temperature.min().values)
        tmin = 5
        tmax = np.ceil(c_ds.temperature.max().values)
        smin = 20
        smax = np.ceil(c_ds.salinity.max().values)

        proj=ccrs.AlbersEqualArea()
        self.plt_proj=ccrs.PlateCarree()
        proj_dict = {"projection": self.plt_proj}
        fig, axs = plt.subplots(4, 2, figsize=(5.5,5.5), subplot_kw=proj_dict)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15,
                            hspace=0.05, wspace=0.05)

        seasons = ['DJF','MAM','JJA','SON']
        for i, season in enumerate(seasons):
            da = c_ds.sel(season=season).load()
            print (da)
            pt = self.render(axs[i,0], da.temperature, crange=(tmin,tmax))
            ps = self.render(axs[i,1], da.salinity, crange=(smin,smax))

        self.set_cbar(axs[-1,0], pt, r'Surface Temperature ($^{\circ}$C)')
        self.set_cbar(axs[-1,1], ps, r'Surface Salinity (-)')

        plt.savefig('FIGS/en4_gridded_surface_climatology.png', dpi=600)

ps = plot_surface()
ps.plot_gridded_surface_en4()
