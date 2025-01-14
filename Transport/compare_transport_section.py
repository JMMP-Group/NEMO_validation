import xarray as xr
import pandas as pd



def preprocess_annual_cruise_data():
    def open_ds(fn, path):
    
        df = pd.read_csv(path + fn, engine="c")
        ds = xr.Dataset.from_dataframe(df)
        
        return ds
    
    path = "../../analysis_eel_data/data/raw/csv_datagridded/"
    data = open_ds("EELCTDandLADCP_3Dfield.csv", path)

    path = "../../analysis_eel_data/data/raw/csv_ctdgrid/"
    pos = open_ds("EELCTDandLADCP_refpos.csv", path)
    time = open_ds("EELCTDandLADCP_refdate.csv", path)
    
    # unstack dist, year and depth
    data = data.set_coords(["Refdist","Year","Depth"])
    data = data.set_xindex(["Refdist","Year","Depth"])
    data = data.unstack("index")
    
    # set distance to index
    pos = pos.set_coords("Refdist")
    pos = pos.swap_dims({"index":"Refdist"})
    pos = pos.drop_vars("index")
    
    # set time to index
    time = time.set_coords("Year")
    time = time.swap_dims({"index":"Year"})
    time = time.drop_vars("index")
    
    ds = xr.merge([data,pos,time])
    
    print (ds)

preprocess_annual_cruise_data()

