from PythonEnvCfg.config import config
cfg = config() # initialise variables in python

import coast
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import cmocean

class plot_gridded_surface(object):
    """ ploting routine for gridded surface temperature and salinity maps """

    def __init__(self):

        # paths
        clim_path = "/surface_maps/surface_state_climatology.nc"
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_dat = cfg.dn_out + clim_path
        en4_nc = "/surface_maps/en4_gridded_surface_climatology.nc"
        m_nc = "/surface_maps/binned_model_surface_climatology.nc"
        self.en4_grid = cfg.dn_out + en4_nc
        self.model_binned = cfg.dn_out + m_nc

    def set_cbar(self, fig, ax, p, txt):
        """ add colour bar to bottom of row """

        pos = ax.get_position()
        cbar_ax = fig.add_axes([pos.x0, 0.08, 
                                pos.x1 - pos.x0, 0.02])
        cbar = fig.colorbar(p, cax=cbar_ax, orientation='horizontal')
        cbar.ax.text(0.5, -3.0, txt, fontsize=8,
                     rotation=0, transform=cbar.ax.transAxes,
                     va='top', ha='center')

    def render(self, ax, da, crange, diff=False):
        """ render subplot panel with data """

        # select colour map
        if diff:
            cmap = cmocean.cm.balance
        elif da.name == 'temperature':
            cmap = cmocean.cm.thermal
        elif da.name == 'salinity':
            cmap = cmocean.cm.haline

        # render
        p = ax.pcolor(da.longitude_bins, da.latitude_bins, da,
                      vmin=crange[0], vmax=crange[1], cmap=cmap,
                      shading='nearest', transform=self.plt_proj)

        # mask with landmass
        ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
        return p

    def plot_native_model_surface_state_climatology(self, scalar="temperature"):
        """ plot surface state climatology over seasonal means """

        # get data
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
            pt = self.render(axs[i,0], da.temperature, crange=(tmin,tmax))
            ps = self.render(axs[i,1], da.salinity, crange=(smin,smax))

        self.set_cbar(axs[-1,0], pt, r'Surface Temperature ($^{\circ}$C)')
        self.set_cbar(axs[-1,1], ps, r'Surface Salinity (-)')

        plt.savefig('FIGS/surface_climatology.png', dpi=600)

    def plot_gridded_surface_en4(self):
        """ plot en4 data on regular lat lon grid """

        # get data
        c_ds = xr.open_dataset(self.en4_grid)
        c_ds = c_ds.rename(dict(potential_temperature="temperature",
                                practical_salinity="salinity"))

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
            pt = self.render(axs[i,0], da.temperature, crange=(tmin,tmax))
            ps = self.render(axs[i,1], da.salinity, crange=(smin,smax))

        self.set_cbar(axs[-1,0], pt, r'Surface Temperature ($^{\circ}$C)')
        self.set_cbar(axs[-1,1], ps, r'Surface Salinity (-)')

        plt.savefig('FIGS/en4_gridded_surface_climatology.png', dpi=600)

    def validate_surface_climatology(self, scalar="temperature"):
        """
        Compare surface climatology of binned model data against EN4.
        """

        # set labels and colour ranges
        if scalar == "temperature": 
            label = r'Surface Temperature ($^{\circ}$C)'
            d_label = r'Surface Temperature Error ($^{\circ}$C)'
            crange = (5, 25)
            diff_range = (-5, 5)
        if scalar == "salinity": 
            label = r'Surface Salinity ($1 \times 10^3$)'
            d_label = r'Surface Salinity Error ($1 \times 10^3$)'
            crange = (20, 35)
            diff_range = (-5, 5)

        # prepare figure
        proj=ccrs.AlbersEqualArea()
        self.plt_proj=ccrs.PlateCarree()
        proj_dict = {"projection": self.plt_proj}
        fig, axs = plt.subplots(4, 3, figsize=(6.5,5.5), subplot_kw=proj_dict)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15,
                            hspace=0.05, wspace=0.05)

        # get en4 data
        en4 = xr.open_dataset(self.en4_grid)
        en4 = en4.rename(dict(potential_temperature="temperature",
                              practical_salinity="salinity"))

        # get model data
        model = xr.open_dataset(self.model_binned)
        model = model.rename(dict(longitude="longitude_bins",
                                  latitude="latitude_bins"))

        # align dims
        model, en4 = xr.align(model, en4)

        # render each panel
        seasons = ['DJF','MAM','JJA','SON']
        for i, season in enumerate(seasons):
            da_en4 = en4.sel(season=season).load()
            da_mod = model.sel(season=season).load()

            self.render(axs[i,1], da_mod[scalar].T, crange=crange)
            pt = self.render(axs[i,0], da_en4[scalar], crange=crange)
            d_pt = self.render(axs[i,2], da_mod[scalar].T - da_en4[scalar],
                        crange=diff_range, diff=True)

        # add colour bars
        self.set_cbar(fig, axs[-1,0], pt, label)
        self.set_cbar(fig, axs[-1,1], pt, label)
        self.set_cbar(fig, axs[-1,2], d_pt, d_label)

        # save
        plt.savefig("FIGS/surface_" + scalar + "_climatology_validation.png",
                    dpi=600)

    def validate_model_surface_climatology(self):
        """ plot all climatology for temperature and salinity """

        self.validate_surface_climatology("temperature")
        self.validate_surface_climatology("salinity")

class plot_pointwise_surface(object):
    """
    ploting routine for pointwise surface temperature and salinity maps
    """

    def __init__(self):

        # paths
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_path = cfg.dn_out + "surface_maps/"
        self.fn_comp_path = cfg.comp_case["proc_data"] + "surface_maps/"
        en4_nc = "/surface_maps/en4_gridded_surface_climatology.nc"
        m_nc = "/surface_maps/binned_model_surface_climatology.nc"
        self.en4_grid = cfg.dn_out + en4_nc
        self.model_binned = cfg.dn_out + m_nc

    def plot_validate_model_surface_by_season(self, var="temperature",
                                              comp_model=None):
        """
        plot modelled surface temperature and salinity against EN4

        parameters
        ----------
        var: scalar variable to plot
        comp_model: comparison model, if required
        """

        proj=ccrs.PlateCarree()
        self.plt_proj=ccrs.PlateCarree()
        self.proj_dict = {"projection": self.plt_proj}

        if comp_model:
            self.render_validate_model_surface_by_season_two_model(var)
        else:
            self.render_validate_model_surface_by_season_one_model(var)

    def render_validate_model_surface_by_season_one_model(self, var):
        """ plot suface comparison accross two model configurations """

        # initiate plots
        fig, axs = plt.subplots(2,2, figsize=(6.5,6.5), subplot_kw=proj_dict)
        plt.subplots_adjust()

        # get data
        ds_path = self.fn_dat + "near_surface_EN4_bias_by_season.nc"
        da = xr.open_dataset(ds_path)["diff_" + var]

        vmin = -1
        vmax = 1

        def render(ax, da):
            p = ax.scatter(da.longitude, da.latitude, c=da, s=0.2,
                           transform=self.plt_proj, cmap=cmocean.cm.balance,
                           vmin=vmin, vmax=vmax)
            return p

        for i, (season, da) in enumerate(da.groupby("season")):
            ax = axs.flatten()[i]
            p = render(ax,da)
            ax.text(0.5,1.01,season,transform=ax.transAxes)

        for ax in axs.flatten():
            ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        pos0 = axs[1,0].get_position()
        pos1 = axs[1,1].get_position()
        cbar_ax = fig.add_axes([pos0.x0, 0.12, 
                                pos1.x1 - pos1.x0, 0.02])
        cbar = fig.colorbar(p, cax=cbar_ax, orientation='horizontal')
        cbar.ax.text(0.5, -2.8, r"Temperature ($^{\circ}$C)", fontsize=8,
                     rotation=0, transform=cbar.ax.transAxes,
                     va='top', ha='center')

        plt.show()

    def render_validate_model_surface_by_season_two_model(self, var):
        """ plot suface comparison accross two model configurations """

        # initiate plots
        fig, axs = plt.subplots(4,3, figsize=(6.5,6.5), 
                                subplot_kw=self.proj_dict)
        plt.subplots_adjust()

        # get data
        ds_path = self.fn_path + "near_surface_EN4_bias_by_season.nc"
        da = xr.open_dataset(ds_path)["diff_" + var]

        ds_path = self.fn_comp_path + "near_surface_EN4_bias_by_season.nc"
        comp_da = xr.open_dataset(ds_path)["diff_" + var]

        da = da.set_index(id_dim=["time","latitude","longitude"])
        comp_da = comp_da.set_index(id_dim=["time","latitude","longitude"])
        print(da)
 
        diff = da - comp_da

        vmin = -1
        vmax = 1

        def render(ax, da):
            p = ax.scatter(da.longitude, da.latitude, c=da, s=0.2,
                           transform=self.plt_proj, cmap=cmocean.cm.balance,
                           vmin=vmin, vmax=vmax)
            return p

        # Render model
        for i, (season, da_season) in enumerate(da.groupby("season")):
            ax = axs[i,0]
            p = render(ax, da_season)
            ax.text(0.5,1.01, season, transform=ax.transAxes)

        # Render comparison model    
        for i, (season, comp_da_season) in enumerate(comp_da.groupby("season")):
            ax = axs[i,1]
            p = render(ax, comp_da_season)
            ax.text(0.5,1.01, season, transform=ax.transAxes)

        # Render difference
        for i, (season, diff_season) in enumerate(diff.groupby("season")):
            ax = axs[i,2]
            p = render(ax, diff_season)
            ax.text(0.5,1.01, season, transform=ax.transAxes)

        for ax in axs.flatten():
            ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        pos0 = axs[3,1].get_position()
        pos1 = axs[3,2].get_position()
        cbar_ax = fig.add_axes([pos0.x0, 0.12, 
                                pos1.x1 - pos1.x0, 0.02])
        cbar = fig.colorbar(p, cax=cbar_ax, orientation='horizontal')
        cbar.ax.text(0.5, -2.8, r"Temperature ($^{\circ}$C)", fontsize=8,
                     rotation=0, transform=cbar.ax.transAxes,
                     va='top', ha='center')
        plt.show()

ps = plot_gridded_surface()
ps =  plot_pointwise_surface()
ps.plot_validate_model_surface_by_season("temperature", comp_model=True)
ps.plot_validate_model_surface_by_season("salinity", comp_model=True)
