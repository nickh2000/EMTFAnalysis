import numpy as np
from cmath import tan
import math
from ROOT import *
from array import *
import Helper as h
import array


#This is the angle between strips in the CSC's
#Most of the time, we are working in eighth strips (10 / 75 / 8 = 1/60 of a degree)
nominal_pitch = 10. / 75.



#This file sources the z-coordinate of each chamber in order to accurately draw chamber-origins in space
global z_source
z_source  = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/EMTFNtuple_geometry_run2.root"

#This folder contains the LUTs that we want to correct
#These should be the same LUTs used to generate the station_station_phi data
global folder 
folder = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/L1Trigger/L1TMuon/data/emtf_luts/ph_lut_GOOD/"

#This folder contains the output LUTs
#This folder should contain all the same LUTs as "folder", however some of files will be later be changed when you run this script
global output
output = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/L1Trigger/L1TMuon/data/emtf_luts/ph_lut_custom/"

#Obtain this file by running station_station_phi.py on EMTFNTuples generated with the LUTs used in folder (which we intend to correct)
#This file will contain the distribution of station-station delta-phis for every station-transition and every sector
global station_station_phi_data
station_station_phi_data = "/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/alignment/station_station_phi_custom_chamber.root"

#Used to obtain the z-values for the chambers
evt_tree  = TChain('EMTFNtuple/tree')
evt_tree.Add(z_source)
#only need one event
evt_tree.GetEvent()



#The following information is used to extract chamber information from its cscid
#CSCID is required to get a chamber's global phi
MAX_ENDCAP = 2
MAX_STATION = 4 
MAX_RING = 4 
MAX_CHAMBER = 36
MAX_LAYER = 6

MIN_ENDCAP = 1
MIN_STATION = 1 
MIN_RING = 1 
MIN_CHAMBER = 1 
MIN_LAYER = 1


BITS_ENDCAP = 3 
BITS_STATION = 3 
BITS_RING = 3 
BITS_CHAMBER = 6 
BITS_LAYER = 3

MASK_ENDCAP = 7 
MASK_STATION = 7 
MASK_RING = 7,
MASK_CHAMBER = 0x077
MASK_LAYER = 7

START_CHAMBER = BITS_LAYER
START_RING = START_CHAMBER + BITS_CHAMBER
START_STATION = START_RING + BITS_RING
START_ENDCAP = START_STATION + BITS_STATION

LOWER_THETA = 8.5
UPPER_THETA = 45
theta_scale = float((UPPER_THETA - LOWER_THETA) / 128.)



#Get cscdetid from chamber information
#offline ring is in form: ME1b, ME12, ME13, ME1a
def getCSCDetID(iendcap, istation, iring, ichamber, ilayer = 0):
    #cscid ring in form ME1a, ME1b, ME12, ME13
    if station == 1:
        iring = (iring + 1) % 4
        if iring == 0: iring = 4

    return (ilayer & MASK_LAYER) | ((ichamber & MASK_CHAMBER) << START_CHAMBER) | ((iring & MASK_RING) << START_RING) | ((istation & MASK_STATION) << START_STATION) | ((iendcap & MASK_ENDCAP) << START_ENDCAP)
    
#get ring from cscdetid
def get_ring(index):
    if (((index  >> START_STATION) & MASK_STATION) == 1):
        iring = (index >> START_RING) & MASK_RING
        i = iring - 1
        if i == 0: i = 4
        return i
    else:
      return ((index >> START_RING) & MASK_RING)

#get chamber number from cscdetid
def get_chamber(index):
    return ((index >> START_CHAMBER) & MASK_CHAMBER)

#get station number from cscdetid
def get_station(index):
    return ((index >> START_STATION) & MASK_STATION)

#get endcap from cscid
def get_endcap(index):
    return ((index >> START_ENDCAP) & MASK_ENDCAP)

#Didn't use this function, but theoretically could shift all chambers by some common (dx, dy) vector
def nudge(dx, dy):
    global nudge_phi
    for iChamb in range(len(evt_tree.geo_cscDetID)):


        #Get chamber coords
        cscid = evt_tree.geo_cscDetID[iChamb]
        z = evt_tree.geo_cscZ[iChamb]
        theta = evt_tree.geo_cscTh[iChamb]
        phi = evt_tree.geo_cscPh[iChamb]

        r = z * math.tan(theta)
        x = r * math.cos(phi)
        y = r * math.Sin(phi)


        #Add shift to geometry
        x += dx
        y += dy

        r = math.sqrt(x*x + y*y)

        new_theta = math.arctan(r / z)
        new_phi = math.arctan(y / x)
        
        #get changes in radial coordinates
        dtheta = new_theta - theta 
        dphi = new_phi - phi


        #save these 
        station = station(cscid)
        endcap = endcap(cscid)
        ring = ring(cscid)
        sector = (chamber - 3) / 6 + 1
        if chamber < 3: sector = 6
        cscid, sb = get_sector_cscid_and_substation(station, ring, cscid)
        if sb == 1: station = 0

        #save this chamber shift in LUT format
        nudge_phi[endcap][sector][station][cscid] = dphi


evt_tree.GetEntry(0)


#get displacement of sector's origin from its hard-coded value (-22 + 15 + 60 * (sector - 1))
#sector origin is defined as teh origin of the lowest-neighbor chamber in station 2
def getSectorPhiDispDeg(endcap, sector):

    sector = 6 if sector == 1 else sector - 1
    return (phi_disp[endcap - 1][sector - 1][2][2] * nominal_pitch/2)

#Get global sector position by adding its displacement to its hard-coded value
def getGlobalSectorPhi(endcap, sector):
    return getSectorPhiDispDeg(endcap, sector) - 22 + 15 + (60 * (sector - 1))

#Get global chamber by adding phi_init (displacement from sector origin) to the sector's phi
def getGlobalChamberPhi(endcap, sector, station, substation, chamber):
    return getGlobalSectorPhi(endcap, sector) + (nominal_pitch/8.) * phi_init[endcap - 1][sector - 1][0 if substation == 1 else station][chamber - 1]

#Add chamber's theta_init (theta displacement from sector origin) to sector theta
def getGlobalTheta(endcap, sector, station, substation, chamber):
    return theta_init[endcap - 1][sector - 1][0 if substation == 1 else station][chamber - 1] * theta_scale + LOWER_THETA

#Convert original LUT to array format in order to draw/correct it
'''These chambers are organized in the following format:
        Endcaps (+1, -1)
            Sectors (1-6)
                "Stations" (Subsector 1, subsector 2, Station 2, Station 3, Station 4)
                    Chambers:
                        chamberID(station 1, subsector 1) - | 15 | 7 | 8 | 9 |
                                                            | 14 | 4 | 5 | 6 |
                                                            | 13 | 1 | 2 | 3 | 
                                                            | 16 | 10| 11| 12|
                                                        - 15, 14, 13, and 16 are neighbor chambers
                                                        - subsector 2 is that same but no neighbors (12 chambers total)
                        RCSCID(station 1, subsector 1)    - | 9 | 7 | 8 | 9 |
                                                            | 6 | 4 | 5 | 6 |
                                                            | 3 | 1 | 2 | 3 |
                                                            | 3 | 1 | 2 | 3 |
                        chamberID (stations 2 - 4) -         | 11 | 4 | 5 | 6 | 7 | 8 | 9 |
                                                        |   10   |   1   |   2   |   3   |
                        
'''
def file_to_array():


    #Phi Init: takes phi difference between chamber origin (wire 0, strip 0/MAX) and sector origin
    phi_init_flat = []
    with open(folder + 'ph_init_neighbor.txt') as f:
        #Read in text file and turn it into an array
        phi_init_flat = f.read().replace('\n', '')
        phi_init_flat = phi_init_flat.split(" ")

    #take out last element (its just a space)
    phi_init_flat = phi_init_flat[:-1]

    #We count the total number of entries to ensure they match-up with the file
    flat_index = 0
    for endcap in range(2):
        for sector in range(6):
            for s in range(5):
                for ch in range(16):
                    #Station 1, subsector 2 has 12 chambers
                    if s > 0 and ch > 11: continue
                    #Station 2-4 have 11 chambers
                    if s > 1 and ch > 10: continue
                    phi_init[endcap][sector][s][ch] = phi_init_flat[ flat_index ]
                    flat_index += 1


    #Number of terms we saved should equal number of elements in the file
    if flat_index == len(phi_init_flat): print("Phi Init Loaded")




    #Gets difference between preset hard phi and geometry phi of first strip in chamber
    #61 chambers per sector, 6 sectors, 2 endcaps
    phi_disp_flat = []
    with open(folder + 'ph_disp_neighbor.txt') as f:
        phi_disp_flat = f.read().replace('\n', '')
        phi_disp_flat = phi_disp_flat.split(" ")

    phi_disp_flat = phi_disp_flat[:-1]

    flat_index = 0
    for endcap in range(2):
        for sector in range(6):
            for s in range(5):
                for ch in range(16):
                    if s > 0 and ch > 11: continue
                    if s > 1 and ch > 10: continue
                    phi_disp[endcap][sector][s][ch] = phi_disp_flat[ flat_index ]
                    flat_index += 1

    if flat_index == len(phi_disp_flat): print("Phi Disp Loaded")



    #Gets difference between preset hard theta and geometry theta of first strip in chamber
    #61 chambers per sector, 6 sectors, 2 endcaps
    theta_init_flat = []
    with open(folder + 'th_init_neighbor.txt') as f:
        theta_init_flat = f.read().replace('\n', '')
        theta_init_flat = theta_init_flat.split(" ")

    theta_init_flat = theta_init_flat[:-1]

    flat_index = 0
    for endcap in range(2):
        for sector in range(6):
            for s in range(5):
                for ch in range(16):
                    if s > 0 and ch > 11: continue
                    if s > 1 and ch > 10: continue
                    theta_init[endcap][sector][s][ch] = theta_init_flat[ flat_index ]
                    flat_index += 1

    if flat_index == len(theta_init_flat): print("theta init Loaded")


    #Gets difference between preset hard theta and geometry theta of first strip in chamber
    #61 chambers per sector, 6 sectors, 2 endcaps
    theta_disp_flat = []
    with open(folder + 'th_disp_neighbor.txt') as f:
        theta_disp_flat = f.read().replace('\n', '')
        theta_disp_flat = theta_disp_flat.split(" ")

    theta_disp_flat = theta_disp_flat[:-1]

    flat_index = 0
    for endcap in range(2):
        for sector in range(6):
            for s in range(5):
                for ch in range(16):
                    if s > 0 and ch > 11: continue
                    if s > 1 and ch > 10: continue
                    theta_disp[endcap][sector][s][ch] = theta_disp_flat[ flat_index ]
                    flat_index += 1

    if flat_index == len(theta_disp_flat): print("theta disp Loaded")






    #Gets difference between wire theta and wire 0
    theta_lut_flat = []
    with open(folder + 'th_lut_neighbor.txt') as f:
        theta_lut_flat = f.read().replace('\n', '')
        theta_lut_flat = theta_lut_flat.split(" ")


    theta_lut_flat = theta_lut_flat[:-1]

    isME11A = False
    flat_index = 0

    #I decided to 1-index for this LUT in order to more efficiently mirror the coniditionals deciding the max-strips and max-wires
    for endcap in range(1, 3):
        for sector in range(1, 7):
            #index by station and then substation this time
            for s in range(1, 5):
                for sb in range(3):
                    for ch in range(1, 17):
                        #substation 0 does not exist for station 1
                        if s == 1 and sb == 0: continue
                        #stations 2-4 only use substation 0
                        if s > 1 and sb > 0: continue
                        is_neighbor = False
                        isME11A = False

                        #Decide whether a chamber is a neighbor or ME1/1a
                        #Consult the map at the start of this function to better understand what's going on here
                        if s == 1:
                            if ch <= 9: rcscid = ch
                            elif ch <= 12:
                                rcscid = ch - 9
                                isME11A = True
                            elif ch == 13: rcscid = 3
                            elif ch == 14: rcscid = 6
                            elif ch == 15: 
                                rcscid = 9
                                isME11A = True
                            elif ch == 16:
                                rcscid = 3
                                is_ME11A = True
                            if ch > 12:
                                is_neighbor = True
                                rsector =  6 if (sector == 1) else sector - 1
                                rsb = 2
                        else:
                            if ch <= 9:
                                rcscid = ch
                            elif ch == 10: rcscid = 3
                            elif ch == 11: rcscid = 6
                            if ch > 9:
                                is_neighbor = True
                                rsector = 6 if (sector == 1) else sector - 1

                        #Set maxWires and maxStrips depending on chamber-type
                        if s == 1:
                            if rcscid <= 3: maxWire = 48
                            elif rcscid <= 6:
                                maxWire = 64
                            else: maxWire = 32
                        elif s == 2:
                            if rcscid <= 3: maxWire = 112
                            else: maxWire = 64
                        else: 
                            if rcscid <= 3: maxWire = 96
                            else: maxWire = 64

                        if s == 1:
                            if isME11A: maxStrip = 48
                            elif rcscid <= 3: maxStrip = 64
                            elif 6 < rcscid and rcscid <= 9: maxStrip = 64
                            else: maxStrip = 80
                        else: maxStrip = 80

                        
                        #Set cuts on number of chambers per subsector
                        if s == 1 and sb == 0: continue
                        if s == 1 and sb == 2 and ch > 12: continue
                        if s != 1 and ch > 11: continue

                        for wire in range(maxWire):
                            station = s
                            if sb == 1: station -= 1
                            #Store theta for every wire in a chamber
                            theta_lut[endcap - 1][sector - 1][station][ch - 1][wire] = theta_lut_flat[ flat_index ]
                            flat_index += 1
                        #account for padding (128 per chamber) for counting the elements
                        flat_index += 128 - maxWire

    if flat_index == len(theta_lut_flat): print("theta lut Loaded")


    theta_corr_lut_flat = []
    with open(folder + 'th_corr_lut_neighbor.txt') as f:
        theta_corr_lut_flat = f.read().replace('\n', '')
        theta_corr_lut_flat = theta_corr_lut_flat.split(" ")

    #This is only ME1/1 (both subsectors 1 and 2)
    theta_corr_lut_flat = theta_corr_lut_flat[:-1]


    #Identical to theta_lut above, but takes wires 1/6, 3/6, 5/6, and every other strip
    flat_index = 0
    for endcap in range(1, 3):
        for sector in range(1, 7):
            for s in range(1, 5):
                for sb in range(3):
                    for ch in range(1, 17):
                        if s == 1 and sb == 0: continue
                        if s > 1 and sb > 0: continue
                        is_neighbor = False
                        isME11A = False

                        if s == 1:
                            if ch <= 9: rcscid = ch
                            elif ch <= 12:
                                rcscid = ch - 9
                                isME11A = True
                            elif ch == 13: rcscid = 3
                            elif ch == 14: rcscid = 6
                            elif ch == 15: 
                                rcscid = 9
                            elif ch == 16:
                                rcscid = 3
                                isME11A = True
                            if ch > 12:
                                is_neighbor = True
                                rsector =  6 if (sector == 1) else sector - 1
                                rsb = 2
                        else:
                            if ch <= 9:
                                rcscid = ch
                            elif ch == 10: rcscid = 3
                            elif ch == 11: rcscid = 9
                            if ch > 9:
                                is_neighbor = True
                                rsector = 6 if (sector == 1) else sector - 1
                        
                        if s == 1:
                            if rcscid <= 3: maxWire = 48
                            elif rcscid <= 6:
                                maxWire = 64
                            else: maxWire = 32
                        elif s == 2:
                            if rcscid <= 3: maxWire = 112
                            else: maxWire = 64
                        else: 
                            if rcscid <= 3: maxWire = 96
                            else: maxWire = 64

                        if s == 1:
                            if isME11A: maxStrip = 48
                            elif rcscid <= 3: maxStrip = 64
                            elif 6 < rcscid and rcscid <= 9: maxStrip = 64
                            else: maxStrip = 80
                        else: maxStrip = 80
                        
                        if s == 1 and sb == 2 and ch > 12: continue
                        if s != 1 and ch > 11: continue

                        if s == 1 and rcscid <= 3 and not isME11A:
                            station  = s
                            if sb == 1: station -= 1
                            #For this LUT, we take wire 1/6, 3/6, 5/6
                            for wire in range(3):
                                #Take every other strip
                                for strip in range(0, maxStrip, 2):
                                    theta_corr_lut[endcap - 1][sector - 1][station][ch - 1][wire * maxStrip // 2 + strip // 2] = theta_corr_lut_flat[ flat_index ]
                                    flat_index += 1
                            flat_index += 128 - 3 * maxStrip // 2

    if flat_index == len(theta_corr_lut_flat): print("theta corr_lut Loaded")
    else: 
        print("Mismatch! Length of File: " + str(len(theta_corr_lut_flat)) + ", Length of array: " + str(flat_index))



#Take station, ring, chamber and return cscid (1-16) and substation (0, 1, 2)
def get_sector_cscid_and_substation(station, ring, chamber):    

    if station == 1 or ring == 2:
        sector_chamber = (chamber - 3) % 6 + 1
    else:
        sector_chamber = (chamber - 2) % 3 + 1

    sb = 0

    #ring order is ME1/b, ME1/2, ME1/3, ME1/1a
    if station == 1:
        sb = 1
        if sector_chamber > 3: 
            sector_chamber -= 3
            sb = 2
        if ring == 2:
            sector_chamber += 3
        if ring == 3:
            sector_chamber += 6
        if ring ==4:
            sector_chamber += 9
    else:
        if ring == 2:
            sector_chamber += 3
    
    return (sector_chamber, sb)
        



#Finds station-station phi displacement per sector,station,endcap
'''Algorithm:
     1. Get station-station-dphi distribution of all sectors and station-station transitions from station_station_phi.py's output root file
     2. Per endcap: find station involved with lowest dphi transition^2 sum over its three transitions
     3. This lowest RMS station is the "reference station", which stays stationary and about which the other stations' sectors align
     4. Find shift for the six sectors in each of the other stations, taking their respective dphi to this "reference station"
'''
def correct_from_data():
    #Take average delta-phi for each endcap, station, and sector transition
    #Sector transitions are (1,2), (1,3),(1,4),(2,3),(2,4)(3,4)
    d_phi_shift = np.zeros(shape = (2, 36), dtype=object)
    f  = TFile(station_station_phi_data)
    for endcap in range(2):
        for chamber in range(36):
            #Get station-station-dphi average per transition/sector
            name = 'd_phi_chamber_%s_%d_%d_%d_c0' % ('(1, 2)', endcap, 0, chamber)
            dphi = f.Get(name)

            gauss_fit = dphi.Fit(formula='gaus', option='Q S')

            dphi_average = gauss_fit.Parameters()[1]
            dphi_average = dphi.GetMean()

            #Due to the non-linear nature of the d-phi average, the shift seems to need to be bigger than the average empirically (20%?)
            d_phi_shift[endcap][chamber] = dphi_average * 1.1

    #Return the array with shifts for each sector
    return d_phi_shift


#Add the shifts from correct_from_data() to each sector in the phi_init array and save our custom ph_init_neighbor.txt file
def write_to_file():
    #"station" is [subsector1, subsector2, station 2, station 3, station 4]
    count = 0
    with open(output + "ph_init_neighbor.txt", 'w') as new_file:
        for endcap in range(2):
            for sector in range(1,7):
                for st in range(5):
                    for ch in range(1, 17):
                        #maximum chambers per subsector
                        if st > 0 and ch > 12: continue
                        if st >1 and ch > 11: continue

                        if st == 0: sb = 1
                        elif st == 1: sb = 2
                        else: sb = 0
                        station = 0 if st == 0 else st - 1

                        #determine if a chamber is a neighbor, if so change its rsector to match its "real" sector
                        rsb = sb
                        rsector = sector
                        rcscid = ch
                        is_neighbor = False
                        if ch <= 9: rcscid = ch
                        elif ch <= 12:
                            rcscid = ch - 9
                            isME11A = True
                        elif ch == 13: rcscid = 3
                        elif ch == 14: rcscid = 6
                        elif ch == 15: 
                            rcscid = 9
                        elif ch == 16:
                            rcscid = 3
                            isME11A = True


                        if ch > 12:
                            is_neighbor = True
                            rsector =  6 if (sector == 1) else sector - 1
                            rsb = 2

                        #get the chamber's phi init from the old file
                        old_phi_init = phi_init[endcap][sector - 1][st][ch - 1]

                        #Get the shift obtained from correct_from_data()
                        #may or may not use rsector here (still testing this)

                        if station == 0 and rcscid < 4:
                            sector_chamber = rcscid + 3 * (rsb - 1)
                            chamber = ((rsector - 1) * 6 + 1 + sector_chamber) % 36 + 1
                            shift = d_phi_shift[int(not endcap)][chamber - 1] / 2.
                            print("Endcap: %d, Station: %d, Sector: %d, Chamber: %d, cscid: %d, sector-chamber: %d, shift: %1.2f" % (not endcap, station + 1, rsector, chamber, ch, rcscid + 3*(rsb - 1), shift))
                        else: shift = 0
                        #Add shift to old LUT and save to new file
                        new_phi_init = int(round(old_phi_init + shift))
                        new_file.write(str(new_phi_init) + " ")

                        #new line every 30 entries                        
                        count += 1
                        if count % 30 == 0: new_file.write("\n")




#Array forms for all 6 LUTs
phi_init = np.zeros(shape=(2, 6, 5, 16))
phi_disp = np.zeros(shape=(2, 6, 5, 16))
theta_init = np.zeros(shape=(2, 6, 5, 16))
theta_disp = np.zeros(shape=(2, 6, 5, 16))
theta_corr_lut = np.zeros(shape=(2, 6, 2, 16, 96))
theta_lut = np.zeros(shape=(2, 6, 5, 16, 112))

#Save files to arrays
file_to_array()

#Find sector shifts from data
d_phi_shift = correct_from_data()


#functions for writing shifts to new LUTs and validating with EMTFNtuple file
print(phi_init[0][5][1][2])
write_to_file()