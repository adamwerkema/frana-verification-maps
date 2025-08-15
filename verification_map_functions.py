'''
AUTHOR: 
Adam Werkema

DATE LAST EDITED: 
8/15/2025

WHAT THIS CODE DOES:
Contains functions that make FRANA verifcation maps.

INPUTS:
Start and end dates/times

Outputs:
Interactive map
'''
###############################################################################
# %% IMPORT PACKAGES
import folium, os, datetime, math, json, requests, base64, time
import pandas as pd
import matplotlib.colors as mcolors
import folium.plugins as plugins
import numpy as np
from netCDF4 import Dataset
from PIL import Image
from folium.features import CustomIcon
from random import random

###############################################################################
# %% Make function for grabbing FRANA values for a given lat/lon
def get_frana_value(lat, lon, frana_data):
    lat_index = int((55 - lat)*100)
    lon_index = int((lon + 130)*100)
    return frana_data[lat_index, lon_index]

###############################################################################
# %% GET ASOS DATA

def get_asos(start, end, verbose):

    WAIT_TIME = 2
    MAX_ATTEMPTS = 10

    if verbose:
        print("Grabbing ASOS data...")

    # Make a list of the stations to use
    station_metadata = pd.read_csv("asos_awos_metadata")
    stations = list(station_metadata.loc[station_metadata["ice_sensor"] == True]["name"])
    stations = [station[1:] for station in stations]

    # Make a list of the UTC days to get ASOS data for
    days = pd.date_range(start = datetime.datetime(int(start.split("_")[0]), int(start.split("_")[1]), int(start.split("_")[2])), 
                            end = datetime.datetime(int(end.split("_")[0]), int(end.split("_")[1]), int(end.split("_")[2])),
                            freq = "D", inclusive = "both")
    day_list = []
    for day in days:
        day_list.append(day.strftime('%Y_%m_%d'))

    # Download ASOS data for the days mentioned
    for day in day_list:
        
        attempt = 0
        while attempt < MAX_ATTEMPTS:

            try:
                # Make sure not all requests go to IEM at the same time
                time.sleep(random() * WAIT_TIME)
                
                # Download all ASOS station data for a day
                all_data = pd.read_csv("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?data=wxcodes&data=ice_accretion_1hr&year1=" + 
                                    day[0:4] + "&month1=" + day[5:7] + "&day1=" + day[8:10] +"&year2=" + day[0:4] + "&month2=" + 
                                    day[5:7] + "&day2=" + day[8:10] + 
                                    "&tz=Etc%2FUTC&format=onlycomma&latlon=no&elev=no&missing=M&trace=T&direct=no&report_type=1&report_type=3&report_type=4", 
                                    low_memory = False)
                
                # Exit the loop if data has been collected
                break

            except:
                
                print("There was an error downloading data for " + day + ". Trying again...")
                
                # Add another attempt    
                attempt += 1

        # If MAX_ATTEMPTS was reached and there is still no data, skip station
        if attempt == MAX_ATTEMPTS:
            print(str(attempt) + " unsuccessful attempts to download data " + "for " + day + ". Skipping station.")
            #return
            continue

        # Remove stations that don't have a Goodrich sensor
        all_data = all_data[all_data["station"].isin(stations)]

        # Replace NaNs with M's
        all_data.fillna("M", inplace = True)

        # Remove dates/times that are outside the bounds
        all_data["datetime"] = pd.to_datetime(all_data["valid"])
        start_datetime = datetime.datetime.strptime(start, "%Y_%m_%d_%H")
        end_datetime = datetime.datetime.strptime(end, "%Y_%m_%d_%H")
        all_data.drop(all_data[all_data["datetime"] <= start_datetime].index, inplace = True)
        all_data.drop(all_data[all_data["datetime"] > end_datetime].index, inplace = True)

        # Combine each day of ASOS data
        if day == day_list[0]:
            sorted_data = all_data
        else:
            sorted_data = pd.concat([sorted_data, all_data], ignore_index = True)

    # Now, sort the data by datetime
    sorted_data.sort_values(by = "datetime", inplace = True)
    sorted_data.reset_index(drop = True, inplace=True)

    if verbose:
        print("Done grabbing ASOS data.")
        print("\n")

    ###############################################################################
    if verbose:
        print("QC'ing ASOS data...")

    # Make a list of dictionaries where all the final ASOS data will go
    dict_list = []

    # Calculate how many hourly metars there should be
    max_hourly_metars = (end_datetime - start_datetime).total_seconds()/3600

    # Now, station by station, run a QC check on whether ice accumulation was legit. Place the QC'ed data
    # in a new dataframe
    for station in set(sorted_data["station"]):
        
        # Get data for a single station
        single_station_data = sorted_data[sorted_data["station"] == station].copy()

        # Get station reset time
        reset_time = station_metadata.loc["K" + station == station_metadata["name"], "reset_minute"].item()

        # A few hardcoded reset times for stations that have changed their reset time:
        if (station == "PRN") & (start_datetime > datetime.datetime(2025,1,1)):
            reset_time = "00"
        if station == "NKT":
            if (start_datetime > datetime.datetime(2025,3,1,19)) & (start_datetime < datetime.datetime(2025,3,5,23)):
                reset_time = "54"
            elif start_datetime > datetime.datetime(2025,3,5,23):
                reset_time = "55"
        if (station == "AQV") & (start_datetime > datetime.datetime(2023,7,7,18)):
            reset_time = "47"
        if (station == "BKB") & (start_datetime > datetime.datetime(2023,12,7,11)):
            reset_time = "49"

        # Remove duplicate rows (for counting hourly totals)
        single_station_data_no_duplicates = single_station_data.drop_duplicates(subset='valid', keep='first')

        # Count the number of hourly metars
        hourly_metar_count = (single_station_data_no_duplicates["valid"].str[14:] == reset_time).sum()

        # If the station doesn't have ice accretion, continue
        if (sum(single_station_data["ice_accretion_1hr"] != "M") == 0) & (sum(single_station_data["wxcodes"].str.contains("FZRA")) == 0) & (sum(single_station_data["wxcodes"].str.contains("FZDZ")) == 0):
            dict_list.append({"station" : "K" + station,
                            "ice" : 0,
                            "minutes_of_fzra" : 0,
                            "hours_without_data" : max_hourly_metars-hourly_metar_count})
            continue

        single_station_data["fzra_min"] = 0
        single_station_data["other_min"] = 0

        # Reset index
        single_station_data.reset_index(inplace=True, drop=True)

        # Run QC check on whether FZRA/FZDZ has occured since the reset time. Also count the number of freezing rain minutes
        for i in range(len(single_station_data)):

            # First check whether this is an hourly time. If not, move on.
            if single_station_data["valid"][i][14:] != reset_time:
                continue

            # Create start time variable
            start_t = single_station_data["datetime"][i] - datetime.timedelta(hours=1)

            # Now, continue to step backward through the observations until we get past the start time.
            # Also make sure that we don't go past the earliest time in the dataset
            go_back_counter = 1
            while (i - go_back_counter >= 0) and (single_station_data["datetime"][i - go_back_counter] >= start_t):
                
                wxcode = single_station_data["wxcodes"][i - go_back_counter]

                # If FZRA/FZDZ occurred
                if ("FZRA" in wxcode) | ("FZDZ" in wxcode):

                    # Calculate duration of FZRA/FZDZ
                    duration = single_station_data["datetime"][i - go_back_counter + 1] - single_station_data["datetime"][i - go_back_counter]
                    single_station_data.loc[i, "fzra_min"] += duration.seconds / 60

                if ("RA" in wxcode) | ("DZ" in wxcode) | ("UP" in wxcode) | ("SN" in wxcode) | ("PL" in wxcode):
                    single_station_data.loc[i, "other_min"] = 1

                go_back_counter += 1

        # Drop rows that are not hourly
        single_station_data.drop(single_station_data[~single_station_data["valid"].str.contains(":" + reset_time)].index, inplace = True)

        # If station does not have data, skip station
        if len(single_station_data) == 0:
            continue
            #print(station + " does not have data. Skipping station.")
        
        # Set ice to 0 if fzra_min is 0
        single_station_data.loc[(single_station_data["fzra_min"] == 0) & (single_station_data["other_min"] == 0), "ice_accretion_1hr"] = 0

        # Convert all "M" to 0 and "T" to 0.0001 in ice
        single_station_data.loc[single_station_data["ice_accretion_1hr"] == "T", "ice_accretion_1hr"] = 0.0001
        single_station_data.loc[single_station_data["ice_accretion_1hr"] == "M", "ice_accretion_1hr"] = 0
        single_station_data["ice_accretion_1hr"] = pd.to_numeric(single_station_data["ice_accretion_1hr"])

        # Place the sums in the summary
        dict_list.append({"station" : "K" + station,
                          "ice" : np.sum(single_station_data["ice_accretion_1hr"]),
                          "minutes_of_fzra" : np.sum(single_station_data["fzra_min"]),
                          "hours_without_data" : max_hourly_metars-hourly_metar_count})
        
    # Place all the dictionaries in a dataframe and save as a CSV
    condensed = pd.DataFrame.from_records(dict_list)
    condensed = condensed.round(decimals = 4)
    condensed.to_csv(end + "_asos", index = False)

    if verbose:
        print("Done QC'ing ASOS data.")
        print("\n")

###############################################################################
# %% GET LSRs
def get_lsr_data(start, end, verbose):

    if verbose:
        print("Grabbing LSR data...")

    # Download reports
    start_str = start[0:4] + "-" + start[5:7] + "-" + start[8:10] + "T" + start[11:13] + ":00Z"
    end_str = end[0:4] + "-" + end[5:7] + "-" + end[8:10] + "T" + end[11:13] + ":00Z"
    
    try:
        data = pd.read_csv("https://mesonet.agron.iastate.edu/cgi-bin/request/gis/lsr.py?wfo[]=ALL&sts=" + start_str + "&ets=" + end_str + "&fmt=csv", on_bad_lines = "skip")

    except:
        data = pd.DataFrame()
        data.to_csv(end + "/lsrs", index = False)
        return

    # Remove any reports that are not freezing rain or ice storm
    data = data[(data["TYPETEXT"] == "FREEZING RAIN") | \
                (data["TYPETEXT"] == "Freezing Rain") | \
                (data["TYPETEXT"] == "ICE STORM") | \
                (data["TYPETEXT"] == "FREEZING DRIZZLE")]
    
    # Remove Alaska and Hawaii reports
    data = data[~data["STATE"].str.contains("AK")]
    data = data[~data["STATE"].str.contains("HI")] 

    # Fix lat/lon rows that have multiple data types    
    data["LAT"] = pd.to_numeric(data["LAT"])
    data["LON"] = pd.to_numeric(data["LON"])

    # Remove miscellaneous rows with missing lat/lon
    data = data[data["LAT"] < 99]

    # Reindex data
    data.reset_index(drop = True, inplace = True)

    # Save to CSV
    data.to_csv(end + "/lsrs", index = False)

    if verbose:
        print("Done getting LSR data.")
        print("\n")

###############################################################################
# %% GET mPING
    
def get_mping(start, end, verbose):

    if verbose:
        print("Grabbing mPING data...")

    # Change the format of the start and end times
    mping_start = start[0:4] + "-" + start[5:7] + "-" + start[8:10] + " " + start[11:13] + ":00:00"
    mping_end = end[0:4] + "-" + end[5:7] + "-" + end[8:10] + " " + end[11:13] + ":00:00"

    # Set up the request and get mPING data
    headers = {'content-type' : 'application/json', 
               'authorization' : 'Token 35a1a5c910636cce026fdaeb96cdcafd3604ce1e'
               }
    payload = {'obtime_gte' : mping_start, 'obtime_lt' : mping_end}
    try:
        response = requests.get('https://mping.ou.edu/mping/api/v2/reports.csv', 
                                headers = headers, 
                                params = payload
                                )
    except:
        data = pd.DataFrame()
        data.to_csv(end + '/mping', index = False)
        return
    
    # Output data to a file
    file = open(end + '/mping','w')
    file.write(response.text)
    file.close()

    if verbose:
        print("Done grabbing mPING data.")
        print("\n")

###############################################################################
# %% MAKE MAP
def make_map(start, end, verbose):

    NO_FRANA = True
    NO_ASOS = True
    NO_LSR = True
    NO_mPING = True

    if verbose:
        print("Grabbing FRANA data...")

    # Make a list of each hour that we need to use for FRANA 1 hr accumulation  
    time_range = pd.date_range(start = datetime.datetime(int(start.split("_")[0]), int(start.split("_")[1]), int(start.split("_")[2]), int(start.split("_")[3])), 
                               end = datetime.datetime(int(end.split("_")[0]), int(end.split("_")[1]), int(end.split("_")[2]), int(end.split("_")[3])),
                               freq = "h", inclusive = "right")

    # Now loop through each hour, grab the FRANA 1 hr accumulation from VMRMS-webtest, and accumulate it
    for time in time_range:
        
        # Parse date and time
        time_str = time.strftime("%Y_%m_%d_%H")
        if verbose:
            print("Getting FRANA 1 hour accumulation for " + time_str + ".")
        year = time_str.split("_")[0]
        month = time_str.split("_")[1]
        day = time_str.split("_")[2]
        hour = time_str.split("_")[3]
        
        # Grab 1-hr FRANA accumulation file from whatever source is available to you.
        attempt = 0
        while attempt < 10:
            try:
                # Grab 1-hr FRANA accumulation file from whatever source is available to you.
                ds = Dataset("insert/path/to/frana/netCDF/file/here")
                break
            except:
                print("Getting FRANA file didn't work for " + time_str + ". Trying again...")
                attempt += 1

        # Skip this hour if FRANA is not present
        if attempt < 10:

            # Open NetCDF file and set negatives to 0
            one_hour_accum = np.array(ds["FRANA_Flat_01h"][:])
            ds.close()
            one_hour_accum[one_hour_accum < 0] = 0

            # Accumulate each file
            if time == time_range[0]:
                total_frana_accum = one_hour_accum.copy()
            else:
                total_frana_accum += one_hour_accum

    if np.max(total_frana_accum) > 0:
        NO_FRANA = False

    if verbose:
        print("Done grabbing FRANA data.")
        print("\n")

    # If there are no FRANA files for this time period, then stop making the map
    if not ('total_frana_accum' in locals()):
        print("No FRANA files for " + end + ". Exiting.")
        return

    # Now begin plotting FRANA data
    if verbose:
        print("Plotting FRANA data on map...")

    # Convert from mm to inches
    total_frana_accum = total_frana_accum / 25.4
    
    # Define the colors and bins we will use to plot on the map
    colors = ['#e0dfe0','#C1C1C1','#A1A1A1','#828282','#E6C7CF','#E1AFB4','#DD969A',\
    '#DA7E7F','#D86564','#AB5671','#9A5371','#885071','#774C71','#664971','#543D7C',\
    '#433973','#35386D','#283666','#1E335E','#294b59','#365c6c','#426c7f','#507e92',\
    '#5e8fa5','#6ca0b8']
    bins = [0.001,0.02,0.03,0.04,0.05,0.08,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,\
    0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,10]

    # Now, loop through each color/bin and add it to the map
    red_array = np.zeros(total_frana_accum.shape)
    green_array = np.zeros(total_frana_accum.shape)
    blue_array = np.zeros(total_frana_accum.shape)
    alpha_array = np.zeros(total_frana_accum.shape)

    for counter in range(len(colors)):
        color = colors[counter]
        bin = bins[counter]
        next_bin = bins[counter + 1]
        red_component, green_component, blue_component = mcolors.hex2color(color)

        indices = (total_frana_accum >= bin) & (total_frana_accum < next_bin)

        red_array[indices] = red_component
        green_array[indices] = green_component
        blue_array[indices] = blue_component
        alpha_array[indices] = 1

    # Now combine the RGBA arrays
    frana_rgba_non_mercator = np.stack([red_array,green_array,blue_array,alpha_array], axis = 2)
    
    # Now, expand the frana rgba array so that we can convert to the mercator projection
    frana_rgba = np.zeros((35000,7000,4))

    # Define the convertion from latitude to cartesian coordinates
    def lat2y(lat):
        return math.log(math.tan(lat*math.pi/180) + 1/math.cos(lat*math.pi/180))

    # Calculate the height of the map in cartesian coordinates
    map_top = lat2y(55)
    map_height = map_top - lat2y(20)

    # Now, loop through the non_mercator array and fill in the mercator array
    current_height_in_pixels = 0
    for i in range(3500):
       
        # Get the current latitude
        current_lat = 55 - i*0.01 - 0.01

        # Calculate how many pixels should be designated to this latitude hundredth
        num_pixels = round(((map_top - lat2y(current_lat)) / map_height) * 35000) - current_height_in_pixels

        # Assign n rows from the non_mercator array to the mercator array. n is num_pixels
        frana_rgba[current_height_in_pixels:current_height_in_pixels+num_pixels,:,:] = frana_rgba_non_mercator[i,:,:]

        # Ascend the map by that many pixels
        current_height_in_pixels += num_pixels

    # Now subsample the image so that it won't be too large. If we don't do this, navigating the
    # interactive map will be very slow
    frana_rgba = frana_rgba[::10]

    # Save the image to a png
    image = Image.fromarray((frana_rgba * 255).astype(np.uint8), 'RGBA')
    image.save(end +"/frana.png")

    # Now initialize the verification map and plot on it
    verification_map = folium.Map(location = [37.5, -95], zoom_start = 4, tiles = None)
    folium.TileLayer("OpenStreetMap", show = True, control = False).add_to(verification_map)
    layer1 = folium.FeatureGroup(name='FRANA')
    folium.raster_layers.ImageOverlay(image = end +"/frana.png", bounds=[[20, -130], [55, -60]], opacity = 0.8).add_to(layer1)
    layer1.add_to(verification_map)

    if verbose:
        print("Done plotting FRANA data.")
        print("\n")

    ###############################################################################
    if verbose:
        print("Plotting ASOS on map...")

    # Read in ASOS data from event
    data = pd.read_csv(end + "_asos")

    # If there is ASOS data for this time period, plot it
    if len(data) != 0:

        # Grab metadata for lat and lons
        metadata = pd.read_csv("asos_awos_metadata")

        # Make lat/lon coodinates for each station
        data["lat"] = 0
        data["lon"] = 0
        for station in data["station"]:
            data.loc[data["station"] == station, "lat"] = float((metadata[metadata["name"] == station]["lat"]).iloc[0])
            data.loc[data["station"] == station, "lon"] = float((metadata[metadata["name"] == station]["lon"]).iloc[0])

        # Make colormap
        marker_bins=[0.00,0.0001,0.02,0.03,0.04,0.05,0.08,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,\
        0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,3]
        marker_cpool=["#008000",'#E0E0E0','#C1C1C1','#A1A1A1','#828282','#E6C7CF','#E1AFB4','#DD969A',\
        '#DA7E7F','#D86564','#AB5671','#9A5371','#885071','#774C71','#664971','#543D7C',\
        '#433973','#35386D','#283666','#1E335E','#294b59','#365c6c','#426c7f','#507e92',\
        '#5e8fa5','#6ca0b8']

        # Add markers to the map with color based on ice accumulation
        layer2 = folium.FeatureGroup(name='ASOS', show = False)
        for index, row in data.iterrows():

            # Get FRANA value. Round to T if necessary
            frana_value = get_frana_value(row["lat"], row["lon"], total_frana_accum)
            if (frana_value > 0) & (frana_value < 0.005):
                frana_value = "Trace"
            elif frana_value == 0:
                frana_value = "None"
            else:
                frana_value = str(round(frana_value,2)) + " in"

            # Place FRANA value back into the dataframe for reference
            data.loc[index, "frana"] = frana_value

            # Get ASOS value. Round to T if necessary
            asos_value = row['ice']
            if (asos_value > 0) & (asos_value < 0.005):
                asos_value = "Trace"
            elif asos_value == 0:
                asos_value = "None"
            else:
                asos_value = str(round(asos_value,2)) + " in"

            # Determine the color of the ASOS station
            bin_index = np.searchsorted(marker_bins, row["ice"], side = "right") - 1
            hex_color = marker_cpool[bin_index]

            # Make markers white if FZRA was reported, but no accumulation
            if (row["minutes_of_fzra"] > 0) and (row["ice"] == 0):
                hex_color = "#FFFFFF"
                asos_value = "Value not provided."
            
            # Calculate difference between FRANA and ASOS
            frana_difference = get_frana_value(row["lat"], row["lon"], total_frana_accum)
            asos_difference = float(row['ice'])
            difference = frana_difference - asos_difference
            if asos_value == "Value not provided.":
                difference = "Unknown"
            elif ((asos_value == "Trace") & (frana_value == "Trace")) | ((asos_value == "None") & (frana_value == "None")):
                difference = "None"
            elif ((asos_value == "Trace") & (frana_value == "None")) | ((asos_value == "None") & (frana_value == "Trace")):
                difference = "Trace"
            elif (difference >= -0.005) & (difference < 0.005):
                difference = "None"
            else:
                difference = str(round(difference,2)) + " in"

            # Make the popup content. Change the color based on overprediction or underprediction
            name_line = row['station']
            duration_line_1 = "FZRA duration: "
            duration_line_2 = str(row['minutes_of_fzra']) + " min"
            ice_sensor_line_1 = "Ice sensor: "
            ice_sensor_line_2 = asos_value
            frana_line_1 = "FRANA ice: "
            frana_line_2 = frana_value
            difference_line_1 = "Difference: "
            difference_line_2 = difference

            if (difference == "Unknown") | (difference == "None"):
                color = "black"
            elif (frana_difference - asos_difference) > 0:
                color = "red"
                if difference != "Trace":
                    difference_line_2 = "+" + difference_line_2
            elif (frana_difference - asos_difference) < 0:
                color = "blue"

            # Create custom HTML string for the popup
            if row["hours_without_data"] == 0:
                popup_html = f"""
                    <div>
                        <strong style="font-size: larger;">{name_line}</strong><br>
                        <strong>{ice_sensor_line_1}</strong> {ice_sensor_line_2}<br>
                        <strong>{frana_line_1}</strong> {frana_line_2}<br>
                        <strong>{difference_line_1}</strong> <span style="color:{color};">{difference_line_2}</span><br>
                        <strong>{duration_line_1}</strong> {duration_line_2}<br>
                    </div>
                """

                popup = folium.Popup(popup_html, max_width = 300)
                icon_properties = {"icon_size": [30,30]}
                folium.Marker(location=[row['lat'], row['lon']],
                            popup = popup,
                            icon = plugins.BeautifyIcon(icon = "none", 
                                                        icon_shape="circle",
                                                        background_color=hex_color,
                                                        popup_anchor=(0,-7),
                                                        **icon_properties)
                            ).add_to(layer2)
            else:
                outage_line_1 = "Possible data outage: "
                outage_line_2 = str(int(row["hours_without_data"])) + " hr"
                popup_html = f"""
                    <div>
                        <strong style="font-size: larger;">{name_line}</strong><br>
                        <strong>{ice_sensor_line_1}</strong> {ice_sensor_line_2}<br>
                        <strong>{frana_line_1}</strong> {frana_line_2}<br>
                        <strong>{difference_line_1}</strong> <span style="color:{color};">{difference_line_2}</span><br>
                        <strong>{duration_line_1}</strong> {duration_line_2}<br>
                        <strong>{outage_line_1}</strong> {outage_line_2}<br>
                    </div>
                """

                popup = folium.Popup(popup_html, max_width = 300)
                icon_properties = {"icon_size": [30,30]}
                folium.Marker(location=[row['lat'], row['lon']],
                            popup = popup,
                            icon = plugins.BeautifyIcon(icon = "bolt",
                                                        text_color="black",            # Make sure the bolt is visible
                                                        border_color="black",
                                                        icon_shape="circle",
                                                        background_color=hex_color,
                                                        popup_anchor=(0,-7),
                                                        extra_classes="fa fa-3x",
                                                        **icon_properties)
                            ).add_to(layer2)

            
        layer2.add_to(verification_map)

        NO_ASOS = False

        # Save the edited ASOS file
        data.to_csv(end + "_asos", index = False)

        if verbose:
            print("Done plotting ASOS data.")
            print("\n")

    else:
        print("No ASOS data for " + end + ". Skipping.") 
        print("\n")

    ###############################################################################
    if verbose:
        print("Plotting LSR data...")

    # Read in LSR data from event
    data = pd.read_csv(end + "/lsrs")

    # If there is LSR data for this time period, plot it
    if len(data) != 0:

        # Make colormap
        marker_bins=[0.00,0.001,0.02,0.03,0.04,0.05,0.08,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,\
        0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,10]
        marker_cpool=["#FFFFFF",'#E0E0E0','#C1C1C1','#A1A1A1','#828282','#E6C7CF','#E1AFB4','#DD969A',\
        '#DA7E7F','#D86564','#AB5671','#9A5371','#885071','#774C71','#664971','#543D7C',\
        '#433973','#35386D','#283666','#1E335E','#294b59','#365c6c','#426c7f','#507e92',\
        '#5e8fa5','#6ca0b8']

        # Add LSRs to the map with color based on ice accumulation
        layer3 = folium.FeatureGroup(name='LSRs', show = False)
        for index, row in data.iterrows():
            
            bin_index = np.searchsorted(marker_bins, row["MAG"], side = "right") - 1
            hex_color = marker_cpool[bin_index]

            # Get FRANA value. Round to T if necessary
            frana_value = get_frana_value(row["LAT"], row["LON"], total_frana_accum)
            if (frana_value > 0) & (frana_value < 0.005):
                frana_value = "Trace"
            elif frana_value == 0:
                frana_value = "None"
            else:
                frana_value = str(round(frana_value,2)) + " in"
            
            name_line = "Local Storm Report"
            time_line_1 = "Report time: "
            time_line_2 = str(row["VALID2"][5:]) + " UTC"
            magnitude_line_1 = "Reported ice: "
            magnitude_line_2 = str(row["MAG"]) + " in"
            frana_line_1 = "FRANA ice: "
            frana_line_2 = frana_value
            source_line_1 = "Source: "
            source_line_2 = str(row["SOURCE"])
            remark_line_1 = "Remarks: "
            remark_line_2 = str(row["REMARK"])
            if remark_line_2 == "nan":
                remark_line_2 = "None"

            # Create custom HTML string for the popup
            popup_html = f"""
                <div>
                    <strong style="font-size: larger;">{name_line}</strong><br>
                    <strong>{time_line_1}</strong> {time_line_2}<br>
                    <strong>{magnitude_line_1}</strong> {magnitude_line_2}<br>
                    <strong>{frana_line_1}</strong> {frana_line_2}<br>
                    <strong>{source_line_1}</strong> {source_line_2}<br>
                    <strong>{remark_line_1}</strong> {remark_line_2}<br>
                </div>
            """
            popup = folium.Popup(popup_html, max_width = 200)
            folium.Marker(location=[row['LAT']+.001, row['LON']],
                        popup = popup,
                        icon=CustomIcon(icon_image="interactive_images/" + hex_color + ".png", icon_size=(20, 20))
                        ).add_to(layer3)
        layer3.add_to(verification_map)

        NO_LSR = False

        if verbose:
            print("Done plotting LSR data.")
            print("\n")

    else:
        if verbose:
            print("No LSR data for " + end + ". Skipping.") 
            print("\n")

    ###############################################################################
    if verbose:
        print("Plotting mPING data...")

    # Read in mPING data
    data = pd.read_csv(end + "/mping")

    # If there is LSR data for this time period, plot it
    if len(data) != 0:

        # Reduce data to just freezing precip types
        data = data.loc[(data["description_id"] == 4) | (data["description_id"] == 6) | (data["description_id"] == 48)]

        # Add mPINGs to the map
        layer4 = folium.FeatureGroup(name='mPING', show = False)
        for index, row in data.iterrows():
            
            # Get FRANA value. Round to T if necessary
            if (row["lat"] > 20) & (row["lat"] < 55) & (row["lon"] > -130) & (row["lon"] < -60):
                frana_value = get_frana_value(row["lat"], row["lon"], total_frana_accum)
                if (frana_value > 0) & (frana_value < 0.005):
                    frana_value = "Trace"
                elif frana_value == 0:
                    frana_value = "None"
                else:
                    frana_value = str(round(frana_value,2)) + " in"
            else:
                frana_value = "Unavailable"
            
            name_line = "mPING"
            time_line_1 = "Report time: "
            time_line_2 = row["obtime"][5:7] + "/" + row["obtime"][8:10] + " " + row["obtime"][11:16] + " UTC"
            type_line_1 = "Reported type: "
            if row["description_id"] == 4:
                type_line_2 = "Freezing rain"
            elif row["description_id"] == 6:
                type_line_2 = "Freezing drizzle"
            else:
                type_line_2 = "Freezing rain & ice pellets"
            frana_line_1 = "FRANA ice: "
            frana_line_2 = frana_value

            # Create custom HTML string for the popup
            popup_html = f"""
                <div>
                    <strong style="font-size: larger;">{name_line}</strong><br>
                    <strong>{time_line_1}</strong> {time_line_2}<br>
                    <strong>{type_line_1}</strong> {type_line_2}<br>
                    <strong>{frana_line_1}</strong> {frana_line_2}<br>
                </div>
            """
            popup = folium.Popup(popup_html, max_width = 250)
            folium.Marker(location=[row['lat'], row['lon']],
                        popup = popup,
                        icon=CustomIcon(icon_image="interactive_images/raindrop.png", icon_size=(20, 20))
                        ).add_to(layer4)
        layer4.add_to(verification_map)

        NO_mPING = False

        if verbose:
            print("Done plotting mPING data.")
            print("\n")

    else:
        print("No mPING data for " + end + ". Skipping.") 
        print("\n")

    ###############################################################################

    if verbose:
        print("Adding legend...")

    with open("interactive_images/ice_colorbar_v2.png", "rb") as lf:
                # open in binary mode, read bytes, encode, decode obtained bytes as utf-8 string
                b64_content = base64.b64encode(lf.read()).decode("utf-8")
                widget = plugins.FloatImage(
                    "data:image/png;base64,{}".format(b64_content),
                    bottom = 0,
                    left=0,
                    width = "750px"
                )
               
    verification_map.add_child(widget)

    if verbose:
        print("Done adding legend.")
        print("\n")

    ###############################################################################
    if verbose:
        print("Adding state borders...")

     # Add state borders
    def style_function(feature):
        return {
            'fillColor': 'transparent',
            'color': 'black',  # Set border color to black
            'weight': 1,       # Border thickness
        }
    folium.GeoJson("interactive_images/state_borders.geojson", 
                   name = "State borders", 
                   style_function = style_function,
                   show = True,
                   control=False,
                   ).add_to(verification_map)
    
    if verbose:
        print("Done adding state borders.")
        print("\n")

    ###############################################################################
    if verbose:
        print("Adding county borders...")

     # Add state borders
    def style_function(feature):
        return {
            'fillColor': 'transparent',
            'color': 'black',  # Set border color to black
            'weight': 0.25,       # Border thickness
        }
    folium.GeoJson("interactive_images/counties.geojson", 
                   name = "County borders", 
                   style_function = style_function,
                   show = False,
                   control=True,
                   ).add_to(verification_map)
    
    if verbose:
        print("Done adding county borders.")
        print("\n")

    ###############################################################################    
    if verbose:
        print("Adding title to map...")

    # Calculate duration between start and end
    date_format = '%Y_%m_%d_%H'
    date1 = datetime.datetime.strptime(start, date_format)
    date2 = datetime.datetime.strptime(end, date_format)
    duration = date2 - date1
    duration_hours = duration.total_seconds() / 3600

    # Make title
    title = str(round(duration_hours)) + "h period ending "
    if end[5] == "0":
        title += end[6]
    else:
        title += end[5:7]
    title += "-"
    if end[8] == "0":
        title += end[9]
    else:
        title += end[8:10]
    title += "-" + end[0:4] + " "
    if end[-2] == "0":
        title += end[-1]
    else:
        title += end[-2:]
    title += " UTC"

    # Make html box with title
    title_html = f'''
                <div style="position: fixed; 
                            top: 10px; left: 50%; transform: translateX(-50%);
                            width: 400px; height: 40px; 
                            background-color: rgba(255, 255, 255, 1); z-index:1000;
                            display: flex; justify-content: center; align-items: center;
                            text-align: center; font-size:20px;
                            border: 2px solid rgba(0, 0, 0, 0.25);
                            border-radius: 5px;"> 
                {title}
                </div>
                '''

    # Add the HTML overlay to the map
    verification_map.get_root().html.add_child(folium.Element(title_html))
    if verbose:
        print("Done adding title to map.")
        print("\n")

    ###############################################################################    
    if verbose:
        print("Adding layer control...")

    # Add layer control and save the map
    folium.LayerControl(collapsed=False).add_to(verification_map)

    if verbose:
        print("Done adding layer control.")
        print("\n")

    ###############################################################################
    if verbose:
        print("Saving map as html file...")

    # Make a location to save maps to (if the location doesn't already exist)
    if not os.path.exists(end[0:4] + "/" + end[5:7] + "/"):
        os.system("mkdir -p " + end[0:4] + "/" + end[5:7] + "/")

    # Save the map as an html file
    verification_map.save(end[0:4] + "/" + end[5:7] + "/" + end + ".html")

    if verbose:
        print("Done saving map as html file.")
        print("\n")

    ###############################################################################
    if verbose:
        print("Cleaning up files...")

    # Remove all the FRANA data from the directory (clean up files we won't use anymore)
    os.system("rm -r " + end)
    
    # Let user know that the map is finished
    print("Done with " + end + ".")
    print("\n")