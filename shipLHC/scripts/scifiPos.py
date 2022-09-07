import ROOT,os
import rootUtils as ut
import shipunit as u
import numpy as np
import time
h={}

ROOT.gROOT.SetBatch(True)
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=False)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=False)
parser.add_argument("--offline", dest="OffMode", default=False, action="store_true")
parser.add_argument("-clusID", dest="ClusterID", required=False, default=0)
parser.add_argument("-procID", dest="ProcID", required=False, default=None, type=int)
options = parser.parse_args()

if options.OffMode:
     xroot_prefix = 'root:://eosuser.cern.ch/'
else:
     xroot_prefix = ''
pathGeofile = xroot_prefix+'/eos/user/a/aiulian/sim_snd/numu_sim_activeemu_withcrisfiles_25_July_2022/1/'
pathSim = xroot_prefix+'/eos/user/a/aiulian/sim_snd/numu_sim_activeemu_withcrisfiles_25_July_2022/'
pathPlots = '/afs/cern.ch/work/f/falicant/public/matching/histo_match/'
pathText = '/afs/cern.ch/work/f/falicant/public/matching/text_match/'
nameGeofile = 'geofile_full.Genie-TGeant4.root'
nameSim = 'sndLHC.Genie-TGeant4_dig.root'
geoFile = pathGeofile+nameGeofile
if options.ProcID is not None:
     simFile = pathSim+str(options.ProcID+1)+'/'+nameSim
else:
     simFile = pathSim+str(1)+'/'+nameSim

#mc_geoFile = '/home/fabio/Simulations_sndlhc/numu_sim_activeemu_withcrisfiles_25_July_2022/1/geofile_full.Genie-TGeant4.root'

start_time = time.time()

fgeo = ROOT.TFile.Open(geoFile)
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
#load geo dictionary
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
# -----Create geometry----------------------------------------------
import shipLHC_conf as sndDet_conf
import SndlhcGeo
geo = SndlhcGeo.GeoInterface(geoFile)
lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])
run = ROOT.FairRunSim()
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()
nav = ROOT.gGeoManager.GetCurrentNavigator()

scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')

zBins = [289, 299, 302, 312, 315, 325, 328, 338, 341, 351, 354, 364]
from array import array
zBins=array('f', zBins)

#from os.path import exists
#cbmsim = ROOT.TChain('cbmsim')
#mc_file_path = xroot_prefix + '/eos/user/a/aiulian/sim_snd/numu_sim_activeemu_withcrisfiles_25_July_2022/'
#mc_file_path = '/home/fabio/Simulations_sndlhc/numu_sim_activeemu_withcrisfiles_25_July_2022/'

#n_files_to_read = 1
#n_files_read = 0
#for i in range (n_files_to_read):
     #file_name = mc_file_path+str(i+1)+'/sndLHC.Genie-TGeant4_dig.root'
     #if not exists(file_name):
     #     continue
     #this_read = cbmsim.Add(file_name)
     #if this_read > 0:
     #     n_files_read += 1
#eventTree = cbmsim
#print('n files read = ', n_files_read)

f=ROOT.TFile.Open(simFile)
eventTree = f.cbmsim

targetRangeX = [-50, 0]
targetRangeY = [10, 60]
targetRangeZ = [270, 370]
wallRangeX = [-46, -6]
wallRangeY = [14, 54]

single = False

if single:
     histoFile = ROOT.TFile('/home/fabio/Simulations_sndlhc/numu_sim_activeemu_withcrisfiles_25_July_2022/histoSciFi_s.root', 'recreate')
elif options.ProcID is not None:
     histoFileName = pathPlots+str(options.ClusterID)+'_histoSciFi_'+str(options.ProcID)+'.root'
     histoFile = ROOT.TFile(histoFileName, 'recreate')
     fitFileName = pathText+str(options.ClusterID)+'_fitSciFi_'+str(options.ProcID)+'.txt'
     statFileName = pathText+str(options.ClusterID)+'_statSciFi_'+str(options.ProcID)+'.txt'
     fitFile = open(fitFileName, "w+")
     statFile = open(statFileName, "w+")
     #fitFile.write('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}  {5:>10}  {6:>10}'.format('event','x_bar','x_rad','y_bar','y_rad','x_hits', 'y_hits')+'\n')
     def addToFile(addPrint):
          fitFile.write(addPrint+'\n')
nBin = 50
sizeBin = 50/nBin
binRange = 3
ut.bookCanvas(h,'digi_dis','digi_hiys distribution',1200, 500, cx=2,cy=1)
ut.bookHist(h,'digi_x','x position ; x [cm]',nBin,targetRangeX[0],targetRangeX[1])
ut.bookHist(h,'digi_y','y position ; y [cm]',nBin,targetRangeY[0],targetRangeY[1])
if single:
     ut.bookCanvas(h,'2d_map','2d map',800, 600, cx=1,cy=1)
     ut.bookHist(h,'digi','2d digi; x[cm];y [cm]',nBin,targetRangeX[0],targetRangeX[1],nBin,targetRangeY[0],targetRangeY[1])
else:
     ut.bookHist(h,'res_x','x residuals; x [cm]',100,-20,20)
     ut.bookHist(h,'res_y','y residuals; y [cm]',100,-20,20)
     ut.bookHist(h,'res_d','distance residuals; [cm]',100,0,25)
print('SciFi start')
A,B = ROOT.TVector3(),ROOT.TVector3()
ccCount=0
Nev = eventTree.GetEntries()
for i_event, event in enumerate(eventTree):
     if Nev<0: break
     if options.OffMode and Nev%100==0:
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
     if single:
          print('event->', i_event)
     ccCount += 1
     h['digi_x'].Reset()
     h['digi_y'].Reset()
     if single:
          h['digi'].Reset()
     for aHit in event.Digi_ScifiHits:
          if not aHit.isValid(): continue
          scifiDet.GetSiPMPosition(aHit.GetDetectorID(),A,B)
          digiStation = int(aHit.GetDetectorID()/1000000)
          if digiStation == motherWall+1:
               if aHit.isVertical():
                    h['digi_x'].Fill(A[0])
                    if single:
                         for bin in np.arange(0, B[1]-A[1],sizeBin):
                              h['digi'].Fill(A[0], B[1]-bin)
               else:
                    h['digi_y'].Fill(B[1])
                    if single:
                         for bin in np.arange(0, B[0]-A[0],sizeBin):
                              h['digi'].Fill(A[0]+bin, B[1])
     mcX = ROOT.TPolyLine()
     mcX.SetPoint(0, nuX, 0)
     mcX.SetPoint(1, nuX, 5)
     mcX.SetLineColor(2)
     mcX.SetLineWidth(2)
     mcY = ROOT.TPolyLine()
     mcY.SetPoint(0, nuY, 0)
     mcY.SetPoint(1, nuY, 5)
     mcY.SetLineColor(2)
     mcY.SetLineWidth(2)
     mc = ROOT.TGraph()
     mc.SetPoint(1, nuX, nuY)
     mc.SetMarkerStyle(29)
     mc.SetMarkerColor(2)
     mc.SetMarkerSize(2)

     h['digi_dis'].cd(1)
     maxBinX = h['digi_x'].GetMaximumBin()
     maxBinY = h['digi_y'].GetMaximumBin()
     maxX = h['digi_x'].GetBinCenter(maxBinX)
     maxY = h['digi_y'].GetBinCenter(maxBinY)
     fitX = h['digi_x'].Fit('gaus','SQ','',maxX-binRange*sizeBin,maxX+binRange*sizeBin)
     fitStatusX = int(fitX)
     fitY = h['digi_y'].Fit('gaus','SQ','',maxY-binRange*sizeBin,maxY+binRange*sizeBin)
     fitStatusY = int(fitY)
     fitMeanX = fitX.Parameter(1)
     fitVarX = fitX.Parameter(2)
     fitMeanY = fitY.Parameter(1)
     fitVarY = fitY.Parameter(2)
     quantile = array('d',[0.5])
     if fitStatusX!=0 or fitMeanX < wallRangeX[0] or fitMeanX > wallRangeX[1]:
          medianX = array('d',[0])
          #print('x fit not converged')
          h['digi_x'].GetQuantiles(1, medianX, quantile)
          fitMeanX = medianX[0]
     if fitStatusY!=0 or fitMeanY < wallRangeY[0] or fitMeanY > wallRangeY[1]:
          medianY = array('d',[0])
          #print('y fit not converged')
          h['digi_y'].GetQuantiles(1, medianY, quantile)
          fitMeanY = medianY[0]
     dist = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
     h['digi_x'].Draw()
     mcX.Draw('same')
     h['digi_dis'].cd(2)
     h['digi_y'].Draw()
     mcY.Draw('same')
     h['digi_dis'].Update()
     if single:
          h['2d_map'].cd(1)
          h['digi'].Draw('COLZ')
          mc.Draw('sameP')
          h['2d_map'].Update()
          histoFile.cd()
          h['digi_dis'].Write()
          histoFile.cd()
          h['2d_map'].Write()

     #stats = h['digi_x_'+str(wall+1)].FindObject('stats')
     #stats.SetOptFit(111)
     if not single:
          h['res_x'].Fill(fitMeanX-nuX)
          h['res_y'].Fill(fitMeanY-nuY)
          h['res_d'].Fill(dist)
     if single:
          print('nu vertex ({:.6f}, {:.6f})'.format(nuX, nuY))
          print('hit - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
          rc = input("hit return for next event or q for quit: ")
          if rc=='q': break
     #print('event {} wall {} res ({:.6f}, {:.6f} -> {:.6f})'.format(i_event, motherWall, fitMeanX-nuX, fitMeanY-nuY, dist))
     if not single:
          addToFile('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}  {5:>10}  {6:>10}'.format(i_event, round(fitMeanX, 2), round(fitVarX, 2), round(fitMeanY, 2), round(fitVarY, 2),h['digi_x'].GetEntries(),h['digi_y'].GetEntries()))
print("Total number of event: {}".format(eventTree.GetEntries()))
print("Selected number of event: {}".format(ccCount))
statFile.write("Total number of event: {}\n".format(eventTree.GetEntries()))
statFile.write("Selected number of event: {}\n".format(ccCount))

if not single:
     ut.bookCanvas(h,'res','residuals',1000, 2500, cx=1,cy=3)
     h['res'].SetCanvasSize(1000, 1500)
     h['res'].SetWindowSize(1050, 1200)
     #h['res'].cd(1)
     #fitResX = h['res_x'].Fit('gaus','SQ','',-10,10)
     #fitResY = h['res_y'].Fit('gaus','SQ','',-10,10)
     #resMeanX = fitResX.Parameter(1)
     #resVarX = fitResX.Parameter(2)
     #resMeanY = fitResY.Parameter(1)
     #resVarY = fitResY.Parameter(2)
     #print('resMeanX {}, resVarX {}, resMeanY {}, resVarY {}'.format(resMeanX, resVarX, resMeanY, resVarY))
     #statFile.write('resMeanX {}, resVarX {}, resMeanY {}, resVarY {}\n'.format(resMeanX, resVarX, resMeanY, resVarY))
     #eff0x = 0
     #eff0y = 0
     #eff1x = h['res_x'].GetEntries()
     #eff1y = h['res_x'].GetEntries()
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
     #statFile.write('efficiency x: {} {} {}\n'.format(effx, eff0x, eff1x))
     #statFile.write('efficiency y: {} {} {}\n'.format(effy, eff0y, eff1y))
     #h['res_x'].Draw()
     #h['res'].cd(2)
     #h['res_y'].Draw()
     #h['res'].cd(3)
     #h['res_d'].Draw()
     #h['res'].Print(mc_file_path+'scifi_res.png')
     histoFile.cd()
     h['res_x'].Write()
     h['res_y'].Write()
     h['res_d'].Write()

histoFile.Write()
histoFile.Close()
fitFile.close()
statFile.close()
print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')
