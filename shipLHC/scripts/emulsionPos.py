import ROOT,os
import rootUtils as ut
import shipunit as u
import numpy as np
import ctypes


h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=False)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
options = parser.parse_args()

fgeo = ROOT.TFile.Open(options.path+options.geoFile)
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
#load geo dictionary
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
# -----Create geometry----------------------------------------------
import shipLHC_conf as sndDet_conf
import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geoFile)
lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])

run = ROOT.FairRunSim()
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()
nav = ROOT.gGeoManager.GetCurrentNavigator()

scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
emulsionDet = ROOT.EmulsionDet()

zBins = [289, 299, 302, 312, 315, 325, 328, 338, 341, 351, 354, 364]
from array import array
zBins=array('f', zBins)

#if options.runNumber>0: 
#    f=ROOT.TFile.Open(options.path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
#    eventTree = f.rawConv
#else:
#    f=ROOT.TFile.Open(options.fname)
#    eventTree = f.cbmsim

from os.path import exists

cbmsim = ROOT.TChain('cbmsim')
mc_file_path = '/home/fabio/Simulations/numu_sim_25_July_2022/'
n_files_to_read = 10
n_files_read = 0
for i in range (n_files_to_read):
     file_name = mc_file_path+str(i)+'/inECC_sndLHC.Genie-TGeant4_dig.root'
     if not exists(file_name):
          continue
     this_read = cbmsim.Add(file_name)
     if this_read > 0:
          n_files_read += 1
eventTree = cbmsim

targetRangeX = [-50, 0]
targetRangeY = [10, 60]
targetRangeZ = [270, 370]

#fitFile = open("/home/fabio/Simulations/numu_sim_25_July_2022/eml_fits.txt", "w+")
#fitFile.write('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}'.format('event','x_bar','x_rad','y_bar','y_rad')+'\n')

#def addToFile(addPrint):
#    fitFile.write(addPrint+'\n')
    
#def sigPos(Nev = -1):
Nev = -1
nBins = 50
sizeBin = 50/nBins
binRange = 3
ut.bookCanvas(h,'eml_map','emulsion map',1200, 1200, cx=1,cy=1)
ut.bookHist(h, 'eml_xy', '2d emulsion map', nBins, targetRangeX[0], targetRangeX[1], nBins, targetRangeY[0], targetRangeY[1])
#ut.bookHist(h, 'eml_z', 'z map', 100, targetRangeZ[0], targetRangeZ[1])
ut.bookCanvas(h,'res','residuals',1000, 2500, cx=1,cy=3)
h['res'].SetCanvasSize(1000, 1500)
h['res'].SetWindowSize(1050, 1200)
ut.bookHist(h,'res_x','x residuals; x [cm]',50,-40,40)
ut.bookHist(h,'res_y','y residuals; y [cm]',50,-40,40)
ut.bookHist(h,'res_d','distance residuals; [cm]',50,0,50)
A, B = ROOT.TVector3(), ROOT.TVector3()
ccCount=0
srX = srY = srD = 0
if Nev < 0: Nev = eventTree.GetEntries()
for event in eventTree:
     if Nev<0: break
     Nevent = eventTree.GetEntries()-Nev
     if Nev%100==0:
               print('events to go:', Nev)
     Nev=Nev-1
     motherWall = -1
     cc = 0
     incoming_nu = event.MCTrack[0]
     outgoing_l = event.MCTrack[1]
     if incoming_nu.GetPdgCode()==14 and outgoing_l.GetPdgCode()==13:
          cc = 1
          mX = incoming_nu.GetStartX()
          mY = incoming_nu.GetStartY()
          mZ = incoming_nu.GetStartZ()
          ROOT.gGeoManager.FindNode(mX, mY, mZ)
          node = ROOT.gGeoManager.GetPath()
          for wall in range(5):
               if 'Wall_'+str(wall) in node:
                         motherWall = wall
                         scifiplane = '/cave_1/Detector_0/volTarget_1/ScifiVolume{}_{}000000'.format(wall+1, wall+1)
                         nav.cd(scifiplane)
                         currNav = nav.GetCurrentNode()
                         S = currNav.GetVolume().GetShape()
                         ox,oy,oz = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
                         P = array('d',[ox,oy,oz])
                         M = array('d',[0,0,0])
                         nav.LocalToMaster(P,M)
                         nuZ = M[2]
                         pX = incoming_nu.GetPx()
                         pY = incoming_nu.GetPy()
                         pZ = incoming_nu.GetPz()
                         #nuE = incoming_nu.GetEnergy()
                         nuX = (pX/pZ)*(nuZ-mZ)+mX
                         nuY = (pY/pZ)*(nuZ-mZ)+mY
     if motherWall == -1: continue
     if not cc: continue
     emulsion = 0
     h['eml_xy'].Reset()
     #h['eml_z'].Reset()
     for eHit in event.EmulsionDetPoint:
          emulsionID = eHit.GetDetectorID()
          nWall = ctypes.c_int(0)
          nRaw = ctypes.c_int(0)
          nColumn = ctypes.c_int(0)
          nPlate = ctypes.c_int(0)
          emlX = eHit.GetX()
          emlY = eHit.GetY()
          #emlZ = eHit.GetZ()
          emulsionDet.DecodeBrickID(emulsionID, nWall, nRaw, nColumn, nPlate)
          nWall = nWall.value
          nRaw = nRaw.value
          nColumn = nColumn.value
          nPlate = nPlate.value
          #print('wall %d, raw %d, column %d, plate %d' %(int(nWall.value), int(nRaw.value), int(nColumn.value), int(nPlate.value)))
          if nWall == motherWall and nPlate == 60:
               trackID = eHit.GetTrackID()
               if trackID < 0: continue
               if event.MCTrack[trackID].GetMotherId()==0:
                    #print('event -> ', Nevent)
                    emulsion = 1
                    #h['eml_z'].Fill(emlZ)
                    h['eml_xy'].Fill(emlX, emlY)
     if not emulsion: continue     
     ccCount+=1
     #print('event {} wall {} res ({:.6f}, {:.6f} -> {:.6f})'.format(Nevent, motherWall, fitMeanX-nuX, fitMeanY-nuY, dist))

     h['eml_map'].cd(1)
     h['eml_xy'].Draw('COLZ')
     #h['eml_map'].cd(2)
     #h['eml_z'].Draw()
     h['eml_map'].Update()
     barX = h['eml_xy'].GetMean(1)
     barY = h['eml_xy'].GetMean(2)
     radX = h['eml_xy'].GetStdDev(1)
     radY = h['eml_xy'].GetStdDev(2)
     #addToFile('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}'.format(Nevent, round(barX, 2), round(radX, 2), round(barY, 2), round(radY, 2)))
     #errMeanX = h['eml_xy'].GetMeanError(1)
     #errMeanY = h['eml_xy'].GetMeanError(2)

     dist = ROOT.TMath.Sqrt((barX-nuX)**2+(barY-nuY)**2)
     srX += abs((barX-nuX))
     srY += abs((barY-nuY))
     srD += dist
     h['res_x'].Fill(barX-nuX)
     h['res_y'].Fill(barY-nuY)
     h['res_d'].Fill(dist)


     #rc = input("hit return for next event or q for quit: ")
     #if rc=='q': break
print("Total number of event: {}".format(eventTree.GetEntries()))
print("Selected number of event: {}".format(ccCount))
print('srX', srX)
print('srY', srY)
print('srD', srD)

h['res'].cd(1)
fitResX = h['res_x'].Fit('gaus','SQ','',-20,20)
fitResY = h['res_y'].Fit('gaus','SQ','',-20,20)
resMeanX = fitResX.Parameter(1)
resVarX = fitResX.Parameter(2)
resMeanY = fitResY.Parameter(1)
resVarY = fitResY.Parameter(2)
#print(resMeanX, resVarX, resMeanY, resVarY)
eff0x = 0
eff0y = 0
eff1x = h['res_x'].GetEntries()
eff1y = h['res_x'].GetEntries()
for binX in range(h['res_x'].GetXaxis().GetNbins()):
     if h['res_x'].GetBinCenter(binX) < resMeanX - 3 * resVarX:
          eff0x += h['res_x'].GetBinContent(binX)
     elif h['res_x'].GetBinCenter(binX) > resMeanX + 3 * resVarX:
          eff0x += h['res_x'].GetBinContent(binX)
for binY in range(h['res_y'].GetXaxis().GetNbins()):
     if h['res_y'].GetBinCenter(binY) < resMeanY - 3 * resVarY:
          eff0y += h['res_y'].GetBinContent(binY)
     elif h['res_y'].GetBinCenter(binY) > resMeanY + 3 * resVarY:
          eff0y += h['res_y'].GetBinContent(binY)
effx = 1-eff0x/eff1x
effy = 1-eff0y/eff1y
print('efficiency x:', effx, eff0x, eff1x)
print('efficiency y:', effy, eff0y, eff1y)
h['res_x'].Draw()
h['res'].cd(2)
h['res_y'].Draw()
h['res'].cd(3)
h['res_d'].Draw()
h['res'].Print('/home/fabio/Simulations/numu_sim_25_July_2022/eml_res_1.png')

#fitFile.close()
#sigPos()