import numpy as np
from cmath import tan
import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h
import argparse
import array

nominal_pitch = 10. / 75.

evt_tree  = TChain('EMTFNtuple/tree')
z_source = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/EMTFNtuple_geometry_run2.root"
evt_tree.Add(z_source)
folder = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/L1Trigger/L1TMuon/data/emtf_luts/ph_lut_v3_data/"

evt_tree.GetEvent()

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

#offline ring is in form: ME1b, ME12, ME13, ME1a
def getCSCDetID(iendcap, istation, iring, ichamber, ilayer = 0):
    #cscid ring in form ME1a, ME1b, ME12, ME13
    if station == 1:
        iring = (iring + 1) % 4
        if iring == 0: iring = 4

    return (ilayer & MASK_LAYER) | ((ichamber & MASK_CHAMBER) << START_CHAMBER) | ((iring & MASK_RING) << START_RING) | ((istation & MASK_STATION) << START_STATION) | ((iendcap & MASK_ENDCAP) << START_ENDCAP)
    

def get_ring(index):
    if (((index  >> START_STATION) & MASK_STATION) == 1):
        iring = (index >> START_RING) & MASK_RING
        i = iring - 1
        if i == 0: i = 4
        return i
    else:
      return ((index >> START_RING) & MASK_RING)

def get_chamber(index):
    return ((index >> START_CHAMBER) & MASK_CHAMBER)

def get_station(index):
    return ((index >> START_STATION) & MASK_STATION)

def get_endcap(index):
    return ((index >> START_ENDCAP) & MASK_ENDCAP)



evt_tree.GetEntry(0)

def getSectorPhiDispDeg(endcap, sector):

    sector = 6 if sector == 1 else sector - 1
    return (phi_disp[endcap - 1][sector - 1][2][2] * nominal_pitch/2)

def getGlobalSectorPhi(endcap, sector):
    return getSectorPhiDispDeg(endcap, sector) - 22 + 15 + (60 * (sector - 1))

def getGlobalChamberPhi(endcap, sector, station, substation, chamber):
    return getGlobalSectorPhi(endcap, sector) + (nominal_pitch/8.) * phi_init[endcap - 1][sector - 1][0 if substation == 1 else station][chamber - 1]

def getGlobalTheta(endcap, sector, station, substation, chamber):
    return theta_init[endcap - 1][sector - 1][0 if substation == 1 else station][chamber - 1] * theta_scale + LOWER_THETA

def file_to_array():
    phi_init_flat = []
    with open(folder + 'ph_init_neighbor.txt') as f:
        phi_init_flat = f.read().replace('\n', '')
        phi_init_flat = phi_init_flat.split(" ")

    phi_init_flat = phi_init_flat[:-1]

    flat_index = 0
    for endcap in range(2):
        for sector in range(6):
            for s in range(5):
                for ch in range(16):
                    if s > 0 and ch > 11: continue
                    if s > 1 and ch > 10: continue
                    phi_init[endcap][sector][s][ch] = phi_init_flat[ flat_index ]
                    flat_index += 1


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







    theta_lut_flat = []
    with open(folder + 'th_lut_neighbor.txt') as f:
        theta_lut_flat = f.read().replace('\n', '')
        theta_lut_flat = theta_lut_flat.split(" ")


    theta_lut_flat = theta_lut_flat[:-1]

    isME11A = False
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

                        
                        if s == 1 and sb == 0: continue
                        if s == 1 and sb == 2 and ch > 12: continue
                        if s != 1 and ch > 11: continue

                        for wire in range(maxWire):
                            station = s
                            if sb == 1: station -= 1
                            theta_lut[endcap - 1][sector - 1][station][ch - 1][wire] = theta_lut_flat[ flat_index ]
                            flat_index += 1
                        flat_index += 128 - maxWire


    if flat_index == len(theta_lut_flat): print("theta lut Loaded")


    theta_corr_lut_flat = []


    with open(folder + 'th_corr_lut_neighbor.txt') as f:
        theta_corr_lut_flat = f.read().replace('\n', '')
        theta_corr_lut_flat = theta_corr_lut_flat.split(" ")

    #This is only ME1/1 (both subsectors 1 and 2)
    theta_corr_lut_flat = theta_corr_lut_flat[:-1]

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
                            for wire in range(3):
                                for strip in range(0, maxStrip, 2):
                                    theta_corr_lut[endcap - 1][sector - 1][station][ch - 1][wire * maxStrip // 2 + strip // 2] = theta_corr_lut_flat[ flat_index ]
                                    flat_index += 1
                            flat_index += 128 - 3 * maxStrip // 2


    if flat_index == len(theta_corr_lut_flat): print("theta corr_lut Loaded")
    else: 
        print("Mismatch! Length of File: " + str(len(theta_corr_lut_flat)) + ", Length of array: " + str(flat_index))

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
        




def correct_from_data():
    d_phi_average = np.zeros(shape = (2, 3, 36, 6), dtype=object)
    d_phi_index_map = [(1, 2), (1, 3), (1, 4),(2, 3), (2, 4), (3, 4)]

    f  = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/alignment/station_station_dphi_RunTwo.root")
    for endcap in range(2):
        for id in range(6):
            for sector in range(1, 7):
                name = 'd_phi_chamber_%s_%d_%d' % (str(d_phi_index_map[id]), endcap, ring, chamber)
                dphi = f.Get(name)

                gauss_fit = dphi.Fit(formula='gaus', option='Q S')

                dphi_average = gauss_fit.Parameters()[1]
                dphi_average = dphi.GetMean()

                d_phi_average[endcap][sector - 1][id] = dphi_average

                
    

    #this approach is naive
    dphi_shift = np.zeros(shape = (2, 6, 4))
    for endcap in range(2):
        #Count stations 2 and 3 as the same station
        station_dphi = np.zeros(shape = (4))
        for sector in range(6):
            for id, avg in enumerate(d_phi_average[endcap][sector]):
                station1 = d_phi_index_map[id][0]
                station2 = d_phi_index_map[id][1]

                station_dphi[station1 - 1] += avg**2
                station_dphi[station2 - 1] += avg**2

        ref_station = np.where(station_dphi == min(station_dphi))[0][0] + 1
        #ref_station= 2
        for sector in range(6):
            for station in range(1 , 5):
                if station == ref_station: continue
                if station > ref_station:
                    transition = (ref_station, station)
                    id = d_phi_index_map.index(transition)
                    dphi_shift[endcap][sector][station - 1] =  d_phi_average[endcap][sector][id]
                else: 
                    transition = (station, ref_station)
                    id = d_phi_index_map.index(transition)
                    dphi_shift[endcap][sector][station - 1] = -1 * d_phi_average[endcap][sector][id]

                    
        for sector in range(6):
            for id, name in enumerate(d_phi_index_map):
                station1_shift = dphi_shift[endcap][sector][name[0] - 1]
                station2_shift = dphi_shift[endcap][sector][name[1] - 1]
                new_dphi = d_phi_average[endcap][sector][id] + station2_shift - station1_shift

                print("Endcap: %d, Transition: (%d,%d), Sector: %d, Old_Dphi: %d, New_Dphi: %d" % (endcap, name[0], name[1], sector, d_phi_average[endcap][sector][id], new_dphi))
    for endcap in range(2):
        for station in range(4):
            for sector in range(6):
                print("Endcap: %d, Station: %d, Sector: %d, Shift: %1.2f" % (endcap, station, sector, dphi_shift[endcap][sector][station]))

    
    return dphi_shift



def write_to_file():
    #"station" is [subsector1, subsector2, station 2, station 3, station 4]
    count = 0
    with open("/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/L1Trigger/L1TMuon/data/emtf_luts/ph_lut_custom/ph_init_neighbor.txt", 'w') as new_file:
        for endcap in range(2):
            for sector in range(1,7):
                for station in range(5):
                    for chamber in range(16):
                        if station > 0 and chamber >= 12: continue
                        if station >1 and chamber >= 11: continue
                        st = 0 if station == 0 else station - 1

                        rsector = sector
                        if st == 0:
                            if chamber >= 12:
                                is_neighbor = True
                                rsector =  6 if (sector == 1) else sector - 1
                        else:
                            if chamber >= 9:
                                is_neighbor = True
                                rsector = 6 if (sector == 1) else sector - 1

                        old_phi_init = phi_init[endcap][sector - 1][station][chamber]

                        #shift is in eighth-strip, disp is in double_strip
                        shift = d_phi_shift[endcap][sector - 1][st]
                        new_phi_init = old_phi_init + shift


                        new_file.write(str(int(round(new_phi_init))) + " ")
                        count += 1
                        if count % 30 == 0: new_file.write("\n")



phi_init = np.zeros(shape=(2, 6, 5, 16))
phi_disp = np.zeros(shape=(2, 6, 5, 16))
theta_init = np.zeros(shape=(2, 6, 5, 16))
theta_disp = np.zeros(shape=(2, 6, 5, 16))
theta_corr_lut = np.zeros(shape=(2, 6, 2, 16, 96))
theta_lut = np.zeros(shape=(2, 6, 5, 16, 112))

file_to_array()

d_phi_shift = correct_from_data()

chamberPosition_old = np.zeros(shape = (2, 4), dtype=object)
chamberPosition_new = np.zeros(shape = (2, 4), dtype=object)


labels = np.empty((2, 4), dtype = list)
for i in np.ndindex(labels.shape): labels[i] = []

chamberPosition_oldPoint = np.empty((2, 4, 2), dtype=list)
for i in np.ndindex(chamberPosition_oldPoint.shape): chamberPosition_oldPoint[i] = []

chamberPosition_newPoint = np.empty((2, 4, 2), dtype=list)
for i in np.ndindex(chamberPosition_newPoint.shape): chamberPosition_newPoint[i] = []

print(chamberPosition_newPoint[0][0][0])


for endcap in range(2):
    for station in range(4):
        
        chamberPosition_old[endcap][station] = TH2D("chamber_position_old_e%d_s%d" % (endcap, station), '', 1000, -5000, 5000, 1000, -5000, 5000)
        chamberPosition_new[endcap][station] = TH2D("chamber_position_new_e%d_s%d" % (endcap, station), '', 1000, -5000, 5000, 1000, 5000, 5000)
        # chamberPosition_old[endcap][station] = TGraph(144)
        # chamberPosition_new[endcap][station] = TGraph(144)

for iChamb in range(len(evt_tree.geo_cscDetID)):
        
        cscid = evt_tree.geo_cscDetID[iChamb]
        z = evt_tree.geo_cscZ[iChamb]

        station = evt_tree.geo_cscStation[iChamb]
        endcap = evt_tree.geo_cscEndcap[iChamb]
        ring = int(evt_tree.geo_cscRing[iChamb])
        chamber = evt_tree.geo_cscChamber[iChamb]


        if station == 1 or ring == 2:
            sector = int((chamber - 3)) / 6 + 1        
            if chamber < 3: sector = 6
        else:
            sector = int(chamber - 2) / 3 + 1
            if chamber < 2: sector = 6

        
        cscid, sb = get_sector_cscid_and_substation(station, ring, chamber)
        
        
        

        theta = getGlobalTheta(endcap, sector, station, sb, cscid)
        original_phi = getGlobalChamberPhi(endcap, sector, station, sb, cscid)
        theta_rad = theta * 3.14159/180.
        original_phi_rad = original_phi * 3.14159/180.

        
        # theta_rad = evt_tree.geo_cscTh[iChamb]
        # original_phi_rad = evt_tree.geo_cscPh[iChamb]
        

        

        r = z * math.tan(theta_rad)
        x = r * math.cos(original_phi_rad)
        old_x = x
        chamberPosition_oldPoint[endcap - 1][station - 1][0].append(x)
        y = r * math.sin(original_phi_rad)
        chamberPosition_oldPoint[endcap - 1][station - 1][1].append(y)


        # if chamber == 2 and station == 2 and ring == 1:
        #     print("Old Geometry---- x: %d, y:%d, endcap: %d, station: %d, ring: %d, chamber: %d, cscid: %d, sector: %d, substation: %d, theta: %d, phi: %d" % (x, y, endcap, station, ring, chamber, cscid, sector, sb, theta, original_phi))
        #     print(getGlobalChamberPhi(endcap, sector, station, sb, cscid))
        #     print(endcap)
        #     print(sector)
        #     print(getGlobalSectorPhi(endcap, sector))
        #     print((nominal_pitch/8.) * phi_init[endcap - 1][sector - 1][0 if sb == 1 else station][cscid - 1])
        sring = str(ring)
        if station == 1:
            if ring == 4: sring = "1a"
            elif ring == 1: sring = "1b"
            else: sring = str(ring)
            
        
        label = "ME%s%d/%s/%d" % ('+' if endcap == 2 else '-', station, str(sring), chamber)
        labels[endcap - 1][station - 1].append(label)
        


        shift = d_phi_shift[endcap - 1][sector - 1][station - 1]

        new_phi = original_phi + shift * nominal_pitch

        new_phi_rad = new_phi * 3.14159 / 180.

        r = float(z) * math.tan(theta_rad)
        x = r * math.cos(new_phi_rad)
        # if shift:
        #     print(old_x)
        #     print(x)
        chamberPosition_newPoint[endcap - 1][station - 1][0].append(x)
        y = r * math.sin(new_phi_rad)
        chamberPosition_newPoint[endcap - 1][station - 1][1].append(y)

        # if x % 600 < 200 or y % 600 < 200:
        #     print("Old Geometry---- shift: %d, x: %d, y:%d, ring: %d, chamber: %d, cscid: %d, sector: %d" % (shift, x, y, ring, chamber, cscid, sector))
        #     print()
        chamberPosition_new[endcap - 1][station - 1].Fill(x, y)

def validate_LUT():
    evt_tree_new = TChain('EMTFNtuple/tree')
    new_coords = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/EMTFNtuple_geometry_custom.root"
    evt_tree_new.Add(new_coords)
    evt_tree_new.GetEntry(0)
    for iChamb in range(len(evt_tree.geo_cscDetID)):
        
        cscid = evt_tree.geo_cscDetID[iChamb]
        z = evt_tree.geo_cscZ[iChamb]

        station = evt_tree.geo_cscStation[iChamb]
        endcap = evt_tree.geo_cscEndcap[iChamb]
        ring = int(evt_tree.geo_cscRing[iChamb])
        chamber = evt_tree.geo_cscChamber[iChamb]

        shift = d_phi_shift[endcap - 1][sector - 1][station - 1] * nominal_pitch / 8.

        phi = evt_tree.geo_cscPh[iChamb]
        sim_phi = phi + shift
        new_phi = evt_tree_new.geo_cscPh[iChamb]

        if abs(new_phi - sim_phi) > 1:
            print(new_phi)
            print(sim_phi)
        else:
            print("Good!")

        theta = evt_tree.geo_cscTh[iChamb]
        new_theta = evt_tree_new.geo_cscTh[iChamb]

        if abs(theta - new_theta) > 1:
            print(new_theta)
            print(theta)
        

        

        

        new_phi = original_phi + shift * nominal_pitch / 8.




for endcap in range(2):
    for station in range(4):

        #Create PDF canvas for overlaying ideal and Run2 pt-reso plots
        canvas = TCanvas(chamberPosition_old[endcap][station].GetName() , chamberPosition_old[endcap][station].GetName(), 700,700)
        n = len(chamberPosition_newPoint[endcap][station][0])
        
        new_graph = TGraph(n, array.array('d', chamberPosition_newPoint[endcap][station][0]), array.array('d', chamberPosition_newPoint[endcap][station][1]))
        old_graph = TGraph(n, array.array('d', chamberPosition_oldPoint[endcap][station][0]), array.array('d', chamberPosition_oldPoint[endcap][station][1] ))

        # leg = TLegend(0.63, 0.7, 0.83, 0.78)
        # leg.AddEntry(chamberPosition_old[endcap][station], 'Old Geometry', 'p')
        # leg.AddEntry(chamberPosition_new[endcap][station], 'New Geometry', 'p')

        # leg.SetMargin(0.5)
        # leg.SetFillStyle(0)
        

        


        gStyle.SetOptStat(0)
        gStyle.SetLegendBorderSize(0)
        gStyle.SetLegendTextSize(0.018);

        gPad.SetGridy(1)
        gPad.SetGridx(1)

        #Configure draw style of these plots
        # chamberPosition_old[endcap][station].SetLineWidth(0)
        # chamberPosition_new[endcap][station].SetLineWidth(0)

        new_graph.SetMarkerColor(kRed)
        old_graph.SetMarkerColor(4)

        new_graph.SetMarkerStyle(20)
        old_graph.SetMarkerStyle(20)


        new_graph.SetMarkerSize(.3)
        old_graph.SetMarkerSize(.5)

        new_graph.SetLineWidth(0)
        old_graph.SetLineWidth(0)

        

        
        #Draw plots to canvas
        old_graph.Draw("AP same")
        new_graph.Draw("P same")

        # mg = TMultiGraph()
        # mg.Add(new_graph, 'lp')
        # mg.Add(old_graph, 'cp')
        # mg.Draw('a')

        # #X and Y axis labels
        # leg.Draw("same")
        # #chamberPosition_old[endcap][station].GetXaxis().SetTitle("(P_{T}^{L1} * q^{L1} * q^{reco}  - P_{T}^{reco}) / P_{T}^{reco}");
        # chamberPosition_old[endcap][station].GetYaxis().SetTitle("X position");
        # chamberPosition_old[endcap][station].GetYaxis().SetLabelSize(.03)
        # chamberPosition_old[endcap][station].GetYaxis().SetTitle("Y position")


        #Save canvas as PDF

        points = TLatex()
        points.SetTextSize(.005)
        points.SetTextFont(42)
        points.SetTextAlign(21)
        for i in range(n):
            label = labels[endcap][station][i]
            x = old_graph.GetPointX(i)
            y = old_graph.GetPointY(i)
            points.DrawText(x, y + .2, label)



        canvas.SaveAs("plots/alignment/pdfs/%s.pdf" % (chamberPosition_old[endcap][station].GetName()))
        del canvas

outfile = TFile("plots/alignment/chamber_position.root", 'recreate')
for endcap in range(2):
    for station in range(4):
        chamberPosition_old[endcap][station].Write()
        chamberPosition_new[endcap][station].Write()

print(d_phi_shift[1][5][3])
print(d_phi_shift[1][5][1])
write_to_file()
validate_LUT()
