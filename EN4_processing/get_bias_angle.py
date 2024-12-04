
class bias_angle(object):
    """
    ploting routine for pointwise surface temperature and salinity maps
    """

    def __init__(self, var):

        # paths
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_path = cfg.dn_out + "surface_maps/"
        self.fn_comp_path = cfg.comp_case["proc_data"] + "surface_maps/"
        en4_nc = "/surface_maps/en4_gridded_surface_climatology.nc"
        m_nc = "/surface_maps/binned_model_surface_climatology.nc"
        self.en4_grid = cfg.dn_out + en4_nc
        self.model_binned = cfg.dn_out + m_nc

        self.var = var

    def get_bias_climatology(self, model)
        # TODO rename to "full depth" not "near full depth"
        fn="near_full_depth_EN4_bias_by_season_by_region.nc"

        # get data
        ds_path = self.fn_path + fn
        da = xr.open_dataset(ds_path)["diff_" + var]

        ds_path = self.fn_comp_path + fn
        da_comp = xr.open_dataset(ds_path)["diff_" + var]

    def get_bias_angle_ds(self, var, depth_var=True):

        for j, region in enumerate(da.region_names.values):
            print ("region: ",  region)
            da_r = da.sel(region_names=region)
            da_comp_r = da_comp.sel(region_names=region)

            multiindex = ["season","time","latitude","longitude"]
            da_r = da_r.set_index(id_dim=multiindex).drop_duplicates("id_dim")
            da_comp_r = da_comp_r.set_index(id_dim=multiindex).drop_duplicates("id_dim")
            da_r, da_comp_r = xr.align(da_r, da_comp_r)

            # get angle
            angle = np.arctan(da_r/da_comp_r)

            seasons = ["DJF","MAM","JJA","SON"]
            colours = ["r","g","b","c"]

            for i, season in enumerate(seasons):
                ax = axs[j,i]

                angle_season = angle.sel(season=season)

                if depth_var:
                    print (depth_var)
                    angle_season = angle_season.swap_dims({"z_dim":"depth"})

                    angle_season_depth = []
                    for i, depth in enumerate(angle.depth):
                        angle_season_depth.append(
                                                angle_season.sel(depth=depth))
                    angle_1d = xr.concat(angle_season_depth, dim="id_dim")


                    s = ax.hist(angle_1d, bins=100,
                                range=[-np.pi/2,np.pi/2])
                    ax.axvline(np.pi/4, c="orange")
                    ax.axvline(-np.pi/4, c="orange")


                ax.set_title(f"{season} {region}")

        for ax in axs[:,0]:
            ax.set_ylabel("Freq.")

        for ax in axs[:,1:].flatten():
            ax.set_yticklabels([])

        for ax in axs.flatten():
            ax.set_xticks([-np.pi/2,-np.pi/4,0,np.pi/4,np.pi/2])

        for ax in axs[:-1,:].flatten():
            ax.set_xticklabels([])

        for ax in axs[-1,:]:
            ax.set_xticklabels(["$-\pi/2$",
                                "$-\pi/4$",
                                "0",
                                "$\pi/4$",
                                "$\pi/2$"])
            ax.set_xlabel(f"{var} bias")

        plt.savefig(f"full_depth_{var}_bias_angle.png", dpi=600)

    def get_bias_angle_hist(self):
        """
        calculate 1-d histogram from bias angle data
        """

    def get_bootstrapped_bias_angle_hist(self):
        """
        Single instance bootstrap sample the bias angle and create histogram
        """



ps =  plot_pointwise_surface()
