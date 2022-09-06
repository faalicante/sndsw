import ROOT,os
import rootUtils as ut
import shipunit as u
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

def GetCharge(pdgcode, pdgdb):
     charge = 0
     if pdgdb.GetParticle(pdgcode):
          charge = pdgdb.GetParticle(pdgcode).Charge()
     elif (pdgcode > 1e+8):
          charge = 1
     return charge

pdgdb = ROOT.TDatabasePDG.Instance()

targetRangeX = [-50, 0]
targetRangeY = [10, 60]
targetRangeZ = [270, 370]
wallRangeX = [-46, -6]
wallRangeY = [14, 54]

#fitFile = open("/home/fabio/Simulations/numu_sim_25_July_2022/eml_fits.txt", "w+")
#fitFile.write('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}  {5:>10}'.format('event','x_bar','x_rad','y_bar','y_rad','hits')+'\n')
#
#def addToFile(addPrint):
#    fitFile.write(addPrint+'\n')
    
Nev = -1
nBins = 50
sizeBin = 50/nBins
binRange = 3
ut.bookCanvas(h,'eml_map','emulsion map',800, 600, cx=1,cy=1)
ut.bookHist(h, 'eml_xy', '2d emulsion map;x [cm]; y [cm]', 1000, targetRangeX[0], targetRangeX[1], 1000, targetRangeY[0], targetRangeY[1])
#ut.bookHist(h, 'eml_z', 'z map', 100, targetRangeZ[0], targetRangeZ[1])
#ut.bookCanvas(h,'res','residuals',1000, 2500, cx=1,cy=3)
#h['res'].SetCanvasSize(1000, 1500)
#h['res'].SetWindowSize(1050, 1200)
#ut.bookHist(h,'res_x','x residuals; x [cm]',100,-10,10)
#ut.bookHist(h,'res_y','y residuals; y [cm]',100,-10,10)
#ut.bookHist(h,'res_d','distance residuals; [cm]',100,0,10)
#ut.bookCanvas(h, 'map', '2d map', 1200, 3600, 1, 3)
#h['map'].SetCanvasSize(1500, 3600)
#h['map'].SetWindowSize(1550, 1200)
#for i in range(1):
#     for j in range(200):
#          ut.bookHist(h, 'eml_map_'+str(i)+str(j), 'eml_tracks_'+str(i)+str(j)+';x[cm];y[cm]', nBins, targetRangeX[0], targetRangeX[1], nBins, targetRangeY[0], targetRangeY[1])
#          ut.bookHist(h, 'scifi_map_'+str(i)+str(j), 'scifi_bar_'+str(i)+str(j)+';x[cm];y[cm]', nBins, targetRangeX[0], targetRangeX[1], nBins, targetRangeY[0], targetRangeY[1])

A, B = ROOT.TVector3(), ROOT.TVector3()
ccCount=0
inCount = 0
outCount = 0
srX = srY = srD = 0
if Nev < 0: Nev = eventTree.GetEntries()
for event in eventTree:
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
                         #scifiplane = '/cave_1/Detector_0/volTarget_1/ScifiVolume{}_{}000000'.format(wall+1, wall+1)
                         #nav.cd(scifiplane)
                         #currNav = nav.GetCurrentNode()
                         #S = currNav.GetVolume().GetShape()
                         #ox,oy,oz = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
                         #P = array('d',[ox,oy,oz])
                         #M = array('d',[0,0,0])
                         #nav.LocalToMaster(P,M)
                         plate60 = '/cave_1/Detector_0/volTarget_1/Wall_{}/Row_0/Brick_0/Emulsion_59'.format(motherWall)
                         nav.cd(plate60)
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
     print('event->', Nevent)
     h['eml_xy'].Reset()
     #h['eml_z'].Reset()
     z60 = M[2]
     for track in event.MCTrack:
          trackMother = track.GetMotherId()
          trackPdg = track.GetPdgCode()
          trackCharge = GetCharge(trackPdg, pdgdb)
          if trackCharge == 0: continue
          if trackMother == 0:
               trX = track.GetStartX()
               trY = track.GetStartY()
               trZ = track.GetStartZ()
               px = track.GetPx()
               py = track.GetPy()
               pz = track.GetPz()
               emlX = (px/pz)*(z60-trZ)+trX
               emlY = (py/pz)*(z60-trZ)+trY
               if wallRangeX[0] < emlX < wallRangeX[1] and wallRangeY[0] < emlY < wallRangeY[1]:
                    emulsion = 1
                    inCount+=1
                    h['eml_xy'].Fill(emlX, emlY)
                    blockIndex = int(ccCount/200)
                    #h['eml_map_'+str(blockIndex)+str(ccCount)].Fill(emlX, emlY)
               else:
                    outCount+=1
                    #print('Error: track out of wall', emlX, emlY)
     if not emulsion: continue     
     ccCount+=1
     if ccCount > 199:break
     mc = ROOT.TGraph()
     mc.SetPoint(1, nuX, nuY)
     mc.SetMarkerStyle(29)
     mc.SetMarkerColor(2)
     mc.SetMarkerSize(2)
     wallBox = ROOT.TPolyLine()
     wallBox.SetPoint(0, wallRangeX[0], wallRangeY[0])
     wallBox.SetPoint(1, wallRangeX[1], wallRangeY[0])
     wallBox.SetPoint(2, wallRangeX[1], wallRangeY[1])
     wallBox.SetPoint(3, wallRangeX[0], wallRangeY[1])
     wallBox.SetPoint(4, wallRangeX[0], wallRangeY[0])
     wallBox.SetLineColor(1)
     wallBox.SetLineWidth(2)
     h['eml_map'].cd(1)
     h['eml_xy'].SetMarkerStyle(20)
     h['eml_xy'].SetMarkerSize(1)
     h['eml_xy'].SetMarkerColor(1)
     h['eml_xy'].Draw('same')
     mc.Draw('sameP')
     wallBox.Draw('same')
     #h['eml_map'].cd(2)
     #h['eml_z'].Draw()
     h['eml_map'].Update()
     barX = h['eml_xy'].GetMean(1)
     barY = h['eml_xy'].GetMean(2)
     radX = h['eml_xy'].GetStdDev(1)
     radY = h['eml_xy'].GetStdDev(2)
     #addToFile('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}  {5:>10}'.format(Nevent, round(barX, 2), round(radX, 2), round(barY, 2), round(radY, 2), h['eml_xy'].GetEntries()))
     #errMeanX = h['eml_xy'].GetMeanError(1)
     #errMeanY = h['eml_xy'].GetMeanError(2)

     #dist = ROOT.TMath.Sqrt((barX-nuX)**2+(barY-nuY)**2)
     #srX += abs((barX-nuX))
     #srY += abs((barY-nuY))
     #srD += dist
     #h['res_x'].Fill(barX-nuX)
     #h['res_y'].Fill(barY-nuY)
     #h['res_d'].Fill(dist)
     print('nu vertex ({:.6f}, {:.6f})'.format(nuX, nuY))
     print('track barycenter/radius x/y', barX, radX, barY, radY)
     print('track - nu position = ({:.6f}, {:.6f})'.format(barX-nuX, barY-nuY))
     rc = input("hit return for next event or q for quit: ")
     if rc=='q': break

#scifiCount = 0
#path = '/home/fabio/Simulations/numu_sim_25_July_2022/'
#scifiFile = path+'scifi_fits.txt'
#with open (scifiFile, 'r') as fs:
#     linesScifi = fs.readlines()
#     for lineScifi in linesScifi:
#          lineScifi = lineScifi.replace('\n', '')
#          lineScifi = lineScifi.split()
#          if lineScifi[0] == 'event': continue
#          scifiBarX = float(lineScifi[1])
#          scifiBarY = float(lineScifi[3])
#          scifiIndex = int(scifiCount/200)
#          h['scifi_map_'+str(scifiIndex)+str(scifiCount)].Fill(scifiBarX, scifiBarY)
#          scifiCount+=1
#          if scifiCount > 199: break

#for i in range(1):
#     #h['eml_map_'+str(i)].SetFillColor(4)
#     #h['scifi_map_'+str(i)].SetFillColor(2)
#     h['map'].cd(1)
#     for j in range(10):
#          color = j%10+1
#          h['eml_map_'+str(i)+str(j)].SetTitle('events')
#
#          h['eml_map_'+str(i)+str(j)].SetMarkerStyle(20)
#          h['scifi_map_'+str(i)+str(j)].SetMarkerStyle(22)
#          if color != 10:
#               h['eml_map_'+str(i)+str(j)].SetMarkerColor(color)
#               h['scifi_map_'+str(i)+str(j)].SetMarkerColorAlpha(color, 0.2)
#          else:
#               h['eml_map_'+str(i)+str(j)].SetMarkerColor(13)
#               h['scifi_map_'+str(i)+str(j)].SetMarkerColorAlpha(13, 0.2)
#          h['eml_map_'+str(i)+str(j)].SetMarkerSize(0.8)
#          h['scifi_map_'+str(i)+str(j)].SetMarkerSize(2)
#          h['eml_map_'+str(i)+str(j)].SetStats(False)
#          h['scifi_map_'+str(i)+str(j)].SetStats(False)
#          h['eml_map_'+str(i)+str(j)].Draw('same')
#          h['scifi_map_'+str(i)+str(j)].Draw('same')
#     #h['map'].cd(i+1).BuildLegend(0.75, 0.85, 0.95, 0.95)
#     h['map'].cd(2)
#     for j in range(100):
#          color = j%10+1
#          h['eml_map_'+str(i)+str(j)].SetMarkerStyle(20)
#          h['scifi_map_'+str(i)+str(j)].SetMarkerStyle(22)
#          if color != 10:
#               h['eml_map_'+str(i)+str(j)].SetMarkerColor(color)
#               h['scifi_map_'+str(i)+str(j)].SetMarkerColorAlpha(color, 0.2)
#          else:
#               h['eml_map_'+str(i)+str(j)].SetMarkerColor(13)
#               h['scifi_map_'+str(i)+str(j)].SetMarkerColorAlpha(13, 0.2)
#          h['eml_map_'+str(i)+str(j)].SetMarkerSize(0.8)
#          h['scifi_map_'+str(i)+str(j)].SetMarkerSize(2)
#          h['eml_map_'+str(i)+str(j)].SetStats(False)
#          h['scifi_map_'+str(i)+str(j)].SetStats(False)
#          h['eml_map_'+str(i)+str(j)].Draw('same')
#          h['scifi_map_'+str(i)+str(j)].Draw('same')
#     #h['map'].cd(i+1).BuildLegend(0.75, 0.85, 0.95, 0.95)
#     h['map'].cd(3)
#     for j in range(200):
#          color = j%10+1
#          h['eml_map_'+str(i)+str(j)].SetMarkerStyle(20)
#          h['scifi_map_'+str(i)+str(j)].SetMarkerStyle(22)
#          if color != 10:
#               h['eml_map_'+str(i)+str(j)].SetMarkerColor(color)
#               h['scifi_map_'+str(i)+str(j)].SetMarkerColorAlpha(color, 0.2)
#          else:
#               h['eml_map_'+str(i)+str(j)].SetMarkerColor(13)
#               h['scifi_map_'+str(i)+str(j)].SetMarkerColorAlpha(13, 0.2)
#          h['eml_map_'+str(i)+str(j)].SetMarkerSize(0.8)
#          h['scifi_map_'+str(i)+str(j)].SetMarkerSize(2)
#          h['eml_map_'+str(i)+str(j)].SetStats(False)
#          h['scifi_map_'+str(i)+str(j)].SetStats(False)
#          h['eml_map_'+str(i)+str(j)].Draw('same')
#          h['scifi_map_'+str(i)+str(j)].Draw('same')
#     #h['map'].cd(i+1).BuildLegend(0.75, 0.85, 0.95, 0.95)
#h['map'].Print('/home/fabio/cernbox/softphys/2dmap.png')
#print("Total number of event: {}".format(eventTree.GetEntries()))
#print("Selected number of event: {}".format(ccCount))
#print("Tracks in/out of target: ", inCount, outCount, outCount/(inCount+outCount))
#
#print('srX', srX)
#print('srY', srY)
#print('srD', srD)
#
#h['res'].cd(1)
#fitResX = h['res_x'].Fit('gaus','SQ','',-20,20)
#fitResY = h['res_y'].Fit('gaus','SQ','',-20,20)
#resMeanX = fitResX.Parameter(1)
#resVarX = fitResX.Parameter(2)
#resMeanY = fitResY.Parameter(1)
#resVarY = fitResY.Parameter(2)
#print('resMeanX {}, resVarX {}, resMeanY {}, resVarY {}'.format(resMeanX, resVarX, resMeanY, resVarY))
#eff0x = 0
#eff0y = 0
#eff1x = h['res_x'].GetEntries()
#eff1y = h['res_y'].GetEntries()
#for binX in range(h['res_x'].GetXaxis().GetNbins()):
#     if h['res_x'].GetBinCenter(binX) < resMeanX - 3 * resVarX:
#          eff0x += h['res_x'].GetBinContent(binX)
#     elif h['res_x'].GetBinCenter(binX) > resMeanX + 3 * resVarX:
#          eff0x += h['res_x'].GetBinContent(binX)
#for binY in range(h['res_y'].GetXaxis().GetNbins()):
#     if h['res_y'].GetBinCenter(binY) < resMeanY - 3 * resVarY:
#          eff0y += h['res_y'].GetBinContent(binY)
#     elif h['res_y'].GetBinCenter(binY) > resMeanY + 3 * resVarY:
#          eff0y += h['res_y'].GetBinContent(binY)
#effx = 1-eff0x/eff1x
#effy = 1-eff0y/eff1y
#print('efficiency x:', effx, eff0x, eff1x)
#print('efficiency y:', effy, eff0y, eff1y)
#h['res_x'].Draw()
#h['res'].cd(2)
#h['res_y'].Draw()
#h['res'].cd(3)
#h['res_d'].Draw()
#h['res'].Print('/home/fabio/cernbox/softpys/eml_res.png')
#
#fitFile.close()