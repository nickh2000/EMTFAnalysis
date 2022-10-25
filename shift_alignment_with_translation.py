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


nominal_pitch = 10. / 75.



evt_tree  = TChain('EMTFNtuple/tree')
z_source = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/EMTFNtuple_geometry_run2.root"
evt_tree.Add(z_source)
folder = "/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/L1Trigger/L1TMuon/data/emtf_luts/ph_lut_run3_2022/"

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
MASK_CHAMBER = x077, MASK_LAYER = 7


START_CHAMBER = BITS_LAYER
START_RING = START_CHAMBER + BITS_CHAMBER
START_STATION = START_RING + BITS_RING
START_ENDCAP = START_STATION + BITS_STATION

LOWER_THETA = 8.5
UPPER_THETA = 45
theta_scale = (UPPER_THETA - LOWER_THETA) / 128


#input ring is in form: ME1b, ME12, ME13, ME1a
def getCSCDetID(iendcap, istation, iring, ichamber, ilayer = 0):
    if station == 1:
        iring = (iring + 1) % 4
        if iring == 0: iring = 4

    return (ilayer & MASK_LAYER) | ((ichamber & MASK_CHAMBER) << START_CHAMBER) | ((iring & MASK_RING) << START_RING) | ((istation & MASK_STATION) << START_STATION) | ((iendcap & MASK_ENDCAP) << START_ENDCAP)
    

def ring(index):
    if (((index  >> START_STATION) & MASK_STATION) == 1):
        iring = (index >> START_RING) & MASK_RING
        i = iring - 1
        if i == 0: i = 4
        return i
    else:
      return ((index >> START_RING) & MASK_RING)

def chamber(index):
    return ((index >> START_CHAMBER) & MASK_CHAMBER)

def station(index):
    return ((index >> START_STATION) & MASK_STATION)

def endcap(index):
    return ((index >> START_ENDCAP) & MASK_ENDCAP)



for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)

    for csc in range(evt_tree.cscSegment_size):
        print(evt_tree.cscSegment_globZ[csc])

def getSectorPhiDispDeg(endcap, sector):
    sector = 6 if sector == 1 else sector - 1
    return (phi_disp[endcap - 1][sector - 1][2][2] * nominal_pitch/2)

def getGlobalSectorPhi(endcap, sector):
    return getSectorPhiDispDeg(endcap, sector) - 22 + 15 + (60 * (sector - 1))

def getGlobalChamberPhi(endcap, sector, station, substation, chamber):

    return getGlobalSectorPhi(endcap, sector) + (nominal_pitch/8) * phi_init[endcap - 1][sector - 1][station if substation != 1 else station - 1][chamber - 1] %360

def getGlobalTheta(endcap, sector, station, substation, chamber):
    return theta_init[endcap - 1][sector - 1][station if substation != 1 else station - 1][chamber - 1] * theta_scale + LOWER_THETA

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
    
    sector_chamber = (chamber - 3) % 6 + 1

    sb = 0
    if station == 1:
        if ring == 1: sb = 1
        if sector_chamber > 3: 
            sector_chamber -= 3
            sb = 2
        if ring == 2:
            sector_chamber += 3
        if ring == 3:
            sector_chamber += 6
    else:
        if ring == 2:
            sector_chamber += 3
    
    return (sector_chamber, sb)
        

def nudge(dx, dy):
    global nudge_phi
    for iChamb in range(len(evt_tree.geo_cscDetID)):
        cscid = evt_tree.geo_cscDetID[iChamb]
        z = evt_tree.geo_cscZ[iChamb]
        theta = evt_tree.geo_cscTh[iChamb]
        phi = evt_tree.geo_cscPh[iChamb]

        r = z * math.tan(theta)
        x = r * math.cos(phi)
        y = r * math.Sin(phi)

        x += dx
        y += dy

        r = math.sqrt(x*x + y*y)

        new_theta = math.arctan(r / z)
        new_phi = math.arctan(y / x)
        

        dtheta = new_theta - theta 
        dphi = new_phi - phi

        station = station(cscid)
        endcap = endcap(cscid)
        ring = ring(cscid)
        sector = (chamber - 3) / 6 + 1
        if chamber < 3: sector = 6
        cscid, sb = get_sector_cscid_and_substation(station, ring, cscid)
        if sb == 1: station = 0

        nudge_phi[endcap][sector][station][cscid]


d_phi_index_map = [(1, 2), (1, 3), (1, 4),(2, 3), (2, 4), (3, 4)]
d_phi_average = np.zeros(shape = (2, 6, 7), dtype=object)

def correct_from_dphi():
    f  = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_unpacked")
    for endcap in range(2):
        for id in range(6):
            for sector in range(1, 7):
                name = 'd_phi_EMTF_%s_%d_%d_c%d' % (str(d_phi_index_map[id]), endcap, sector, 1)
                pos_dphi = f.Get(name)

                gauss_fit = pos_dphi.Fit(formula='gaus', option='Q S')

                pos_dphi_average = gauss_fit.Parameters()[1]

                name = 'd_phi_EMTF_%s_%d_%d_c%d' % (str(d_phi_index_map[id]), endcap, sector, 1)
                neg_dphi = f.Get(name)

                gauss_fit = neg_dphi.Fit(formula='gaus', option='Q S')

                neg_dphi_average = gauss_fit.Parameters()[1]

                d_phi_average[endcap][sector - 1][id] = pos_dphi_average + neg_dphi_average
    

    #this approach is naive
    dphi_shift = np.zeros(2, 6, 4)
    for endcap in range(2):
            for sector in range(6):
                #Count stations 2 and 3 as the same station
                station_dphi = np.zeros(shape = (4))

                for id, avg in d_phi_average[endcap][sector]:
                    station1 = d_phi_index_map[id][0]
                    station2 = d_phi_index_map[id][1]

                    station_dphi[station1 - 1] += avg**2
                    station_dphi[station2 - 1] += avg**2

                ref_station = station_dphi.index(min(station_dphi)) + 1

                for station in range(1 , 5):
                    if station == ref_station: continue
                    if station > ref_station:
                        transition = "(%d, %d)" % (ref_station, station)
                        id = d_phi_index_map.index(transition)
                        dphi_shift[endcap][sector - 1][station] = -1 * d_phi_average[endcap][sector][id]
                    else: 
                        transition = "(%d, %d)" % (station, ref_station)
                        id = d_phi_index_map.index(transition)
                        dphi_shift[endcap][sector][station - 1] = d_phi_average[endcap][sector][id]
                
    return dphi_shift
            

phi_init = np.zeros(shape=(2, 6, 5, 16))
phi_disp = np.zeros(shape=(2, 6, 5, 16))
theta_init = np.zeros(shape=(2, 6, 5, 16))
theta_disp = np.zeros(shape=(2, 6, 5, 16))
theta_corr_lut = np.zeros(shape=(2, 6, 2, 16, 96))
theta_lut = np.zeros(shape=(2, 6, 5, 16, 112))
nuge_phi = np.zeros(shape=(2, 6, 5, 16))
nudge_theta = np.zeros(shape=(2, 6, 5, 16))

file_to_array()

d_phi_shift = correct_from_dphi()

#"station" is [subsector1, subsectoe2, station 2, station 3, station 4]
# count = 0
# with open("/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src/L1Trigger/L1TMuon/data/emtf_luts/ph_lut_custom/phi_disp_neighbor.txt") as new_file:
#     for endcap in range(2):
#         for sector in range(6):
#             for station in range(5):
#                 for chamber in range(16):
#                     if station > 0 and chamber >= 12: continue
#                     if station >1 and chamber >= 11: continue
#                     st = 0 if station == 0 else station - 1
#                     new_phi[endcap][sector][station][chamber] = phi_disp[endcap][sector][station][chamber] + d_phi_shift[endcap][sector][st]
#                     new_file.write(str(new_phi[endcap][sector][station][chamber] + " "))
#                     count += 1
#                     if count % 30 == 0: new_file.write("\n")





chamberPosition_old = TH2D("chamber_position_old", '', 100, -500, 500, 100, -500, 500)
chamberPosition_new = TH2D("chamber_position_new", '', 100, -500, 500, 100, -500, 500)

for iChamb in range(len(evt_tree.geo_cscDetID)):
        cscid = evt_tree.geo_cscDetID[iChamb]
        z = evt_tree.geo_cscZ[iChamb]

        station = station(cscid)
        endcap = endcap(cscid)
        ring = ring(cscid)
        sector = (chamber - 3) / 6 + 1        
        if chamber < 3: sector = 6
        cscid, sb = get_sector_cscid_and_substation(station, ring, cscid)

        theta = getGlobalTheta(endcap, station, sector, sb, cscid)
        original_phi = getGlobalChamberPhi(endcap, sector, station, sb, cscid)
        

        r = z * math.tan(theta)
        x = r * math.cos(original_phi)
        y = r * math.Sin(original_phi)

        chamberPosition_old.Fill(x, y)

        shift = d_phi_shift[endcap - 1][sector - 1][0 if sb == 1 else station][cscid]
        phi_init[endcap - 1][sector - 1][0 if sb == 1 else station][cscid] += shift
        new_phi = getGlobalChamberPhi(endcap, sector, station, sb, cscid)

        r = z * math.tan(theta)
        x = r * math.cos(new_phi)
        y = r * math.Sin(new_phi)

        chamberPosition_new.Fill(x, y)




