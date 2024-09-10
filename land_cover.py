from glob import glob
import numpy as np
import geopandas as gpd
from pyhdf.SD import SD, SDC
import h5py
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap
from scipy.interpolate import interpn
from adjustText import adjust_text
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


land_cover_data = [
    {"Land Cover": "Evergreen Needleleaf Forest", "Class": 1, "Color (in Hex)": "#008000"},
    {"Land Cover": "Evergreen Broadleaf Forest", "Class": 2, "Color (in Hex)": "#00FF00"},
    {"Land Cover": "Deciduous Needleleaf Forest", "Class": 3, "Color (in Hex)": "#99CC00"},
    {"Land Cover": "Deciduous Broadleaf Forest", "Class": 4, "Color (in Hex)": "#99FF99"},
    {"Land Cover": "Mixed Forest", "Class": 5, "Color (in Hex)": "#339966"},
    {"Land Cover": "Closed Shrubland", "Class": 6, "Color (in Hex)": "#993366"},
    {"Land Cover": "Open Shrubland", "Class": 7, "Color (in Hex)": "#FFCC99"},
    {"Land Cover": "Woody Savannas", "Class": 8, "Color (in Hex)": "#CCFFCC"},
    {"Land Cover": "Savannas", "Class": 9, "Color (in Hex)": "#FFCC00"},
    {"Land Cover": "Grasslands", "Class": 10, "Color (in Hex)": "#FF9900"},
    {"Land Cover": "Permanent Wetlands", "Class": 11, "Color (in Hex)": "#006699"},
    {"Land Cover": "Croplands", "Class": 12, "Color (in Hex)": "#FFFF00"},
    {"Land Cover": "Urban and Built-Up", "Class": 13, "Color (in Hex)": "#FF0000"},
    {"Land Cover": "Cropland/Natural Vegetation Mosaic", "Class": 14, "Color (in Hex)": "#999966"},
    {"Land Cover": "Snow and Ice", "Class": 15, "Color (in Hex)": "#FFFFFF"},
    {"Land Cover": "Barren or Sparsely Vegetated", "Class": 16, "Color (in Hex)": "#808080"},
    {"Land Cover": "Water", "Class": 17, "Color (in Hex)": "#000080"},
    # {"Land Cover": "Error", "Class": 255, "Color (in Hex)": "#000000"},
]

class_to_color = {item["Class"]: item["Color (in Hex)"] for item in land_cover_data}

unique_classes = sorted(class_to_color.keys())
hex_colors = [class_to_color[cls] for cls in unique_classes]
rgb_colors = [tuple(int(h[i:i+2], 16)/255. for i in (1, 3, 5)) for h in hex_colors]
IGBP_cmap = ListedColormap(rgb_colors)


MISSING_VALUES = 17

class LandCoverPlotter:
    # for the default parameter, see also:
    # https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pdf
    def __init__(self, MCDC_path, ne_path,
                 pixel_size = 463.312716525, pixel_per_tile = 2400, 
                 uly = 10007554.677, ulx = -20015109.354, 
                 earth_radius = 6371007.181, is_h5file = False):
        self.MCDC_path = MCDC_path
        self.cities = gpd.read_file(ne_path)
        self.pixel_size = pixel_size
        self.pixel_per_tile = pixel_per_tile
        self.uly = uly
        self.ulx = ulx
        self.earth_radius = earth_radius
        self.is_h5file = is_h5file
        
        self.sinusoidal_projection = ccrs.Sinusoidal(
            central_longitude=0, 
            globe=ccrs.Globe(semimajor_axis=self.earth_radius, 
                             semiminor_axis=self.earth_radius))
        self.plate_carree_projection = ccrs.PlateCarree()

    def get_file_name(self, ih, iv):
        return glob(f"{self.MCDC_path}h{ih:02d}v{iv:02d}.*")

    def read_land_cover(self, ih, iv):
        file_names = self.get_file_name(ih, iv)
        if len(file_names) == 0:
            data_array = np.ones([2400, 2400]) * MISSING_VALUES
            return data_array
        else:
            file_name = file_names[0]
        
        if self.is_h5file:
            with h5py.File(file_name, "r") as file:
                # Extract the data array
                data_array = file["/MCD12Q1/Data Fields/LC_Type1"][:]
        else:
            file = SD(file_name, SDC.READ)
            data_array = file.select("LC_Type1")[:]
            file.end()
        return data_array.astype(np.uint8)

    def hv2ullr(self, hid, vid):
        ulx = self.ulx + hid * self.pixel_size * self.pixel_per_tile
        uly = self.uly - vid * self.pixel_size * self.pixel_per_tile
        
        lrx = ulx + self.pixel_size * self.pixel_per_tile
        lry = uly - self.pixel_size * self.pixel_per_tile
        return ulx, uly, lrx, lry
    
    def tile_locator(self, lons: np.array, lats: np.array):
        points = self.sinusoidal_projection.transform_points(
            self.plate_carree_projection, lons, lats)
        xs, ys = points[:, 0], points[:, 1]
        hs = np.floor((xs - self.ulx) / (self.pixel_size * self.pixel_per_tile)).astype(int)
        vs = np.floor((self.uly - ys) / (self.pixel_size * self.pixel_per_tile)).astype(int)
        return np.vstack([hs, vs]).T

    def latlon2landcover(self, lons, lats):
        hvs = self.tile_locator(lons, lats)
        hv_uniques = np.unique(hvs, axis=0)
        ans = np.zeros_like(lats, dtype=np.uint8)
        
        points = self.sinusoidal_projection.transform_points(
            self.plate_carree_projection, lons, lats)[:, :2]

        for h, v in hv_uniques:
            ulx, uly, lrx, lry = self.hv2ullr(h, v)
            data = self.read_land_cover(h, v)
            
            idx = (hvs[:, 0]==h) & (hvs[:, 1]==v)
            xi = points[idx, :]
            
            ans[idx] = interpn(
                (np.linspace(ulx, lrx, 2400), np.linspace(uly, lry, 2400)), 
                data.T, xi, method="nearest", bounds_error=True, fill_value=255
            ).astype(np.uint8)
        return ans
    
    def fix_name(self, old2new: dict):
        for old, new in old2new.items():
            self.cities.loc[self.cities['NAME'] == old, 'NAME'] = new
    
    def city_locator(self, city_name: str):
        return self.cities.loc[self.cities['NAME']==city_name, 
                               ["LONGITUDE", "LATITUDE"]].values
        
    def get_cities(self, lon_min, lon_max, lat_min, lat_max):
        return self.cities.cx[lon_min:lon_max, lat_min:lat_max].sort_values("NATSCALE")

    def show_n_largest_city(self, ax, lon_min, lon_max, lat_min, lat_max, n = 10):
        city_in_range = self.get_cities(lon_min, lon_max, lat_min, lat_max)
        texts = []
        for _, row in city_in_range.iloc[-n:].iterrows():
            name = row["NAME"]
            city_lon, city_lat = row["geometry"].x, row["geometry"].y
            # print(name, lon, lat)
            texts.append(ax.text(city_lon, city_lat, name, 
                transform=self.plate_carree_projection, fontsize=6))
            ax.plot(city_lon, city_lat, 'ro', markersize = 1, 
                    transform=self.plate_carree_projection)

        _ = adjust_text(texts, only_move={'points':'y', 'texts':'y'}, 
                        arrowprops=dict(arrowstyle="->", color='r', lw=0.5))

    def show_cities(self, ax, city_names, lon_min, lon_max, lat_min, lat_max):
        texts = []
        city_in_range = self.get_cities(lon_min, lon_max, lat_min, lat_max)
        for city_name in city_names:
            city_lon, city_lat = city_in_range.loc[self.cities['NAME']==city_name, 
                               ["LONGITUDE", "LATITUDE"]].values[0, :]
            texts.append(ax.text(city_lon, city_lat, city_name, 
                transform=self.plate_carree_projection, fontsize=6))
            ax.plot(city_lon, city_lat, 'ro', markersize = 1, 
                    transform=self.plate_carree_projection)

        _ = adjust_text(texts, only_move={'points':'y', 'texts':'y'}, 
                        arrowprops=dict(arrowstyle="->", color='r', lw=0.5))

    def show_lat_lon(self, ax, lon_delta, lat_delta):
        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, 
                        color='black', alpha=1.0, linestyle=(0, (0.5, 0.5)), 
                        linewidth=0.5,
                        crs = self.plate_carree_projection)
        gl.top_labels = False
        gl.right_labels = False
        # gl.xlocator = mticker.FixedLocator([-15, 0, 15])
        # gl.ylocator = mticker.FixedLocator([45, 50, 55])
        gl.xlocator = mticker.MultipleLocator(lon_delta)
        gl.ylocator = mticker.MultipleLocator(lat_delta)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        return gl