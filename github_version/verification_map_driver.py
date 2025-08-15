'''
AUTHOR: 
Adam Werkema

DATE LAST EDITED: 
8/15/2025

WHAT THIS CODE DOES:
This code makes interactive maps that serve as verification for FRANA. This code serves
as the driver for all the other functions used to make these interactive maps. First, verification 
data is collected from ASOS stations, MPING, and LSRs. Then, FRANA is plotted on the map. Lastly,
all the verification data is plotted on the map.

INPUTS:
Start and end dates/times of verification

Outputs:
Interactive map with FRANA verification
'''

###############################################################################
# %% IMPORT PACKAGES
import sys, argparse, os, pytz, math
from datetime import datetime, timedelta
sys.path.insert(0, "/location/where/verification_map_functions/is/stored")
import verification_map_functions

###############################################################################
# %% USER SETTINGS

# OPTION 1: Get arguments from command line
# parser = argparse.ArgumentParser()
# parser.add_argument("-start", help = "Start yyyy_mm_dd_hh", required = True)
# parser.add_argument("-end", help = "End yyyy_mm_dd_hh", required = True)
# args = parser.parse_args()
# start = args.start
# end = args.end

# OPTION 2: Make start and end dates/times in yyyy_mm_dd_hh
start = "2025_03_28_12"
end = "2025_03_31_12"

# OPTION 3: Use the current date & time to make a 24hr map
# now_utc = datetime.now(pytz.utc)
# seconds_since_midnight_utc = (now_utc - now_utc.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
# interval_seconds = 6 * 3600
# rounded_interval_seconds = math.floor(seconds_since_midnight_utc / interval_seconds) * interval_seconds
# end_rounded_time_utc = now_utc.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(seconds=rounded_interval_seconds)
# start_rounded_time_utc = end_rounded_time_utc - timedelta(days=1)
# start = str(start_rounded_time_utc.year) + "_" + \
#         str(start_rounded_time_utc.month).zfill(2) + "_" + \
#         str(start_rounded_time_utc.day).zfill(2) + "_" + \
#         str(start_rounded_time_utc.hour).zfill(2)
# end = str(end_rounded_time_utc.year) + "_" + \
#         str(end_rounded_time_utc.month).zfill(2) + "_" + \
#         str(end_rounded_time_utc.day).zfill(2) + "_" + \
#         str(end_rounded_time_utc.hour).zfill(2)

###############################################################################
# %% GATHER DATA AND MAKE MAP

# Make a location to place temporary files into
if not os.path.exists(end):
    os.mkdir(end)

verbose = False

verification_map_functions.get_asos(start, end, verbose = verbose)
verification_map_functions.get_lsr_data(start, end, verbose = verbose)
verification_map_functions.get_mping(start, end, verbose = verbose)
verification_map_functions.make_map(start, end, verbose = verbose)