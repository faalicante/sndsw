import ROOT,os
import rootUtils as ut
import shipunit as u
import numpy as np
import ctypes

h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=True)
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

def xPos(detID):
     orientation = (detID//100000)%10
     nStation = 2*(detID//1000000-1)+orientation
     mat = (detID%100000)//10000
     X = detID%1000+(detID%10000)//1000*128
     return [nStation,mat,X]   # even numbers are Y (horizontal plane), odd numbers X (vertical plane)

if options.runNumber>0: 
     f=ROOT.TFile.Open(options.path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
     eventTree = f.rawConv
else:
     f=ROOT.TFile.Open(options.fname)
     eventTree = f.cbmsim

#from os.path import exist
#
#f=ROOT.TFile.Open('/eos/user/c/cvilela/SND_MC_June21/neutrino/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/'+str(fileIndex)+'/sndLHC.Genie-TGeant4_digCPP.root')
#cbmsim = ROOT.TChain('cbmsim')
#mc_file_path = '/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/'
#n_files_to_read = 401
#n_files_read = 0
#for i in range (n_files_to_read):
#     file_name = mc_file_path+str(i)+'/sndLHC.Genie-TGeant4_digCPP.root'
#     if not exist(file_name):
#          continue
#     this_read = cbmsim.Add(file_name)
#     if this_read > 0:
#          n_files_read += 1
#eventTree = cbmsim

def slopes(Nev=-1):
     A,B = ROOT.TVector3(),ROOT.TVector3()
     ut.bookHist(h,'slopesX','slope diffs',1000,-1.0,1.0)
     ut.bookHist(h,'slopesY','slope diffs',1000,-1.0,1.0)
     ut.bookHist(h,'clX','cluster size',20,0.5,20.5)
     ut.bookHist(h,'clY','cluster size',20,0.5,20.5)
# assuming cosmics make straight line
     if Nev < 0: Nev = eventTree.GetEntries()
     for event in eventTree:
          if Nev<0: break
          Nev=Nev-1
          clusters = makeClusters(event)
          xHits = {}
          yHits = {}
          for s in range(5):
               xHits[s]=[]
               yHits[s]=[]
          for aCl in clusters:
               #fill A,B with mean position (energy)
               aCl.GetPosition(A,B)
               vertical = int(aCl.GetFirst()/100000)%10==1
               s = int(aCl.GetFirst()/1000000)-1
               if vertical: 
                    xHits[s].append(ROOT.TVector3(A))
                    rc = h['clX'].Fill(aCl.GetN())
               else: 
                    yHits[s].append(ROOT.TVector3(A))
                    rc = h['clY'].Fill(aCl.GetN())
          proj = {'X':xHits,'Y':yHits}
          for p in proj:
               sls = []
               for s1 in range(5):
                    if len(proj[p][s1]) !=1: continue
                    cl1 = proj[p][s1][0]
                    for s2 in range(s1+1,5):
                         if len(proj[p][s2]) !=1: continue
                         cl2 = proj[p][s2][0]
                         dz = abs(cl1[2]-cl2[2])
                         if dz < 5: continue
                         dzRep = 1./dz
                         m = dzRep*(cl2-cl1)
                         sls.append(m)
               for ix1 in range(len(sls)-1):
                    for ix2 in range(ix1+1,len(sls)):
                         if p=="X": rc = h['slopes'+p].Fill(sls[ix2][0]-sls[ix1][0])
                         if p=="Y": rc = h['slopes'+p].Fill(sls[ix2][1]-sls[ix1][1])                 
     ut.bookCanvas(h,'clusters',' ',1024,768,1,2)
     k=1
     for cl in ['clX', 'clY']:
          h['clusters'].cd(k)
          h[cl].SetTitle('size of ' + cl + ' projection')
          h[cl].Draw()
          h['clusters'].Update()
          k+=1
     ut.bookCanvas(h,'slopes',' ',1024,768,1,2)
     k=1
     for slope in ['slopesX', 'slopesY']:
          h['slopes'].cd(k)
          h[slope].GetXaxis().SetRangeUser(-0.2,0.2)
          h[slope].SetTitle(slope+' projection; delta slope [rad]')
          h[slope].Draw()
          h[slope].Fit('gaus','S','',-0.02,0.02)
          h['slopes'].Update()
          stats = h[slope].FindObject('stats')
          stats.SetOptFit(1111)
          stats.Draw()          
          k+=1

def hitMaps(Nev=-1):
     for mat in range(30):
          ut.bookHist(h,'mat_'+str(mat),'hit map / mat',512,-0.5,511.5)
          ut.bookHist(h,'sig_'+str(mat),'signal / mat',150,0.0,150.)
     N=-1
     if Nev < 0 : Nev = eventTree.GetEntries()
     for event in eventTree:
          N+=1
          if N>Nev: break
          for aHit in event.Digi_ScifiHits:
               if not aHit.isValid(): continue
               X =  xPos(aHit.GetDetectorID())
               rc = h['mat_'+str(X[0]*3+X[1])].Fill(X[2])
               rc  = h['sig_'+str(X[0]*3+X[1])].Fill(aHit.GetSignal(0))
     ut.bookCanvas(h,'hitmaps',' ',1024,768,6,5)
     ut.bookCanvas(h,'signal',' ',1024,768,6,5)
     for mat in range(30):
          tc = h['hitmaps'].cd(mat+1)
          A = h['mat_'+str(mat)].GetSumOfWeights()/512.
          if h['mat_'+str(mat)].GetMaximum()>10*A: h['mat_'+str(mat)].SetMaximum(10*A)
          h['mat_'+str(mat)].Draw()
          tc = h['signal'].cd(mat+1)
          h['sig_'+str(mat)].Draw()

def eventTime():
 Tprev = -1
 freq = 160.316E6
 ut.bookHist(h,'Etime','delta event time; dt [s]',100,0.0,1.)
 ut.bookCanvas(h,'E',' ',1024,2*768,1,2)
 eventTree.GetEvent(0)
 t0 =  eventTree.EventHeader.GetEventTime()/160.E6
 eventTree.GetEvent(eventTree.GetEntries()-1)
 tmax = eventTree.EventHeader.GetEventTime()/160.E6
 ut.bookHist(h,'time','elapsed time; t [s]',1000,t0,tmax)
 for event in eventTree:
    T = event.EventHeader.GetEventTime()
    dT = 0
    if Tprev >0: dT = T-Tprev
    Tprev = T
    rc = h['Etime'].Fill(dT/freq)
    rc = h['time'].Fill( (T/freq-t0))
 tc = h['E'].cd(1)
 rc = h['Etime'].Fit('expo')
 tc.Update()
 stats = h['Etime'].FindObject("stats")
 stats.SetOptFit(111)
 tc = h['E'].cd(2)
 h['time'].Draw()

def mergeSignals(hstore):
  ut.bookHist(hstore,'signalAll','signal all mat',150,0.0,150.)
  for mat in range(30):
    hstore['signalAll'].Add(hstore['sig_'+str(mat)])
  hstore['signalAll'].Scale(1./hstore['signalAll'].GetSumOfWeights())

def signalZoom(smax):
  for mat in range(30):
    h['sig_'+str(mat)].GetXaxis().SetRangeUser(0.,smax)
    tc = h['signal'].cd(mat+1)
    tc.Update()

def beamSpot():
    A,B = ROOT.TVector3(),ROOT.TVector3()
    ut.bookHist(h,'bs','beam spot',100,-100.,10.,100,0.,80.)
    for event in eventTree:
        xMean = 0
        yMean = 0
        w=0
        for d in event.Digi_ScifiHits:
            detID = d.GetDetectorID()
            s = int(detID/1000000)
            modules['Scifi'].GetSiPMPosition(detID,A,B)
            vertical = int(detID/100000)%10==1
            if vertical: xMean+=A[0]
            else: yMean+=A[1]
            w+=1
        rc = h['bs'].Fill(xMean/w,yMean/w)

def clusterDist(Nev=-1):
     ut.bookCanvas(h,'cls_z','clusters z',cx=2,cy=5)
     h['cls_z'].SetCanvasSize(1200, 1500)
     h['cls_z'].SetWindowSize(1250, 1200)
     #ut.bookCanvas(h, 'cls_size', 'hit per cluster', 1024, 720, 1, 1)
     #ut.bookHist(h, 'h_size', 'cluster hits', 20, 0.5, 20.5)
     for wall in range(5):
          #ut.bookCanvas(h,'2d_map'+str(wall),'2d map wall '+str(wall),cx=2,cy=5)
          #h['cls_map'+str(wall)].SetCanvasSize(1200, 2300)
          #h['cls_map'+str(wall)].SetWindowSize(1250, 1200)
          h['hit_z_wall_'+str(wall)] = ROOT.TH1D('hit z wall '+str(wall), 'hit z wall '+str(wall)+'; z [cm]', len(zBins)-1, zBins)
          h['cls_z_wall_'+str(wall)] = ROOT.TH1D('cls z wall '+str(wall), 'cls z wall '+str(wall)+'; z [cm]', len(zBins)-1, zBins)
          #for plane in range(1,6):
          #     ut.bookHist(h,'hit_'+str(wall)+'_'+str(plane),'hit wall '+str(wall)+' plane '+str(plane)+'; x [cm]; y [cm]',50,-50,0,50,10,60)
          #     ut.bookHist(h,'cls_'+str(wall)+'_'+str(plane),'cls wall '+str(wall)+' plane '+str(plane)+'; x [cm]; y [cm]',50,-50,0,50,10,60)
     nInt=0
     nIntWall=0
     nIntBorder=0
     nIntScifi=0
     A,B = ROOT.TVector3(),ROOT.TVector3()
     if Nev < 0: Nev = eventTree.GetEntries()
     totcl=0
     tothit=0
     cccount=0
     for event in eventTree:
          if Nev<0: break
          Nev=Nev-1
          N = eventTree.GetEntries()-Nev
          motherWall = -1
          cc = 0
          ncl=0
          nhit=0
          for mcTrack in event.MCTrack:
               if mcTrack.GetMotherId()==-1 and mcTrack.GetPdgCode()==14:
                    nInt+=1
                    mX = mcTrack.GetStartX()
                    mY = mcTrack.GetStartY()
                    mZ = mcTrack.GetStartZ()
                    ROOT.gGeoManager.FindNode(mX, mY, mZ)
                    node = ROOT.gGeoManager.GetPath()
                    #print('interacting node', node)
                    for wall in range(5):
                         if 'Wall_'+str(wall) in node:
                              nIntWall+=1
                              motherWall = wall
                              #print(eventTree.GetEntries()-Nev, motherWall)
                         elif 'Wallborder_'+str(wall) in node:
                              nIntBorder+=1
                         elif 'ScifiVolume'+str(wall+1) in node:
                              nIntScifi+=1
               elif mcTrack.GetMotherId()==0 and mcTrack.GetPdgCode()==13:
                    cc = 1
          if motherWall == -1: continue
          if cc: 
               cccount +=1
               for wall in range(5):
                    if motherWall==wall:
                         clusters = makeClusters(event)
                         print(N, len(clusters))
                         for aCl in clusters:
                              totcl+=1
                              ncl+=1
                              aCl.GetPosition(A,B)
                              #print(A[0], A[1], A[2], B[0], B[1], B[2])
                              h['cls_z_wall_'+str(wall)].Fill(A[2])
                              #h['h_size'].Fill(aCl.GetN())
                              clStation = int(aCl.GetFirst()/1000000)
                              #print(N, clStation)
                              for plane in range(1, 6):
                                   if plane == clStation:
                                        pass
                                        #h['cls_'+str(wall)+'_'+str(plane)].Fill(A[0], B[1])
                                        #h['cls_'+str(wall)+'_'+str(plane)].Fill(B[0], B[1])

                         prevID = -1
                         for aHit in event.ScifiPoint:
                              partID = aHit.GetTrackID()
                              hitStation = aHit.station()
                              nhit+=1
                              tothit+=1
                              if partID == -2: continue
                              if partID!=prevID:
                                   hitZ = aHit.GetZ()
                                   h['hit_z_wall_'+str(wall)].Fill(hitZ)
                                   for plane in range(1, 6):
                                        if plane == hitStation:
                                             hitX = aHit.GetX()
                                             hitY = aHit.GetY()
                                             #hitEloss = aHit.GetEnergyLoss()
                                             #h['hit_'+str(wall)+'_'+str(plane)].Fill(hitX, hitY)
                                             #h['2d_sig'].Fill(hitX, hitY, hitEloss)
                                             prevID = partID
               print('event -> number of cluster/hit',N, ncl, nhit)
     print('number of mu_neutrino interactions:', nInt)
     print('number of mu_neutrino interactions in walls:', nIntWall)
     print('number of mu_neutrino cc interactions:', cccount)
     print('total number of cluster/hit', totcl, tothit)
     #h['cls_size'].cd()
     #h['h_size'].Draw()
     for wall in range(5):
          h['cls_z'].cd(wall*2+1)
          h['hit_z_wall_'+str(wall)].Draw()
          h['cls_z'].cd((wall+1)*2)
          h['cls_z_wall_'+str(wall)].Draw()
          #for plane in range(1, 6):
          #     h['cls_map'+str(wall)].cd(plane*2-1)
          #     h['hit_'+str(wall)+'_'+str(plane)].Draw('COLZ')
          #     h['cls_map'+str(wall)].cd(plane*2)
          #     h['cls_'+str(wall)+'_'+str(plane)].Draw('COLZ')
     
def singleEvent(start=0,wall = 0,save=False):
     nBins = 50
     #ut.bookCanvas(h,'z_dis','clusters z',1200, 1000, cx=1,cy=2)
     #h['hit_z_wall_'+str(wall)] = ROOT.TH1D('hit z wall '+str(wall), 'number of hits wall '+str(wall)+'; z [cm]', len(zBins)-1, zBins)
     #h['cls_z_wall_'+str(wall)] = ROOT.TH1D('cls z wall '+str(wall), 'number of clusters wall '+str(wall)+'; z [cm]', len(zBins)-1, zBins)
     ut.bookCanvas(h,'2d_map'+str(wall),'2d map wall '+str(wall),cx=2,cy=5)
     h['2d_map'+str(wall)].SetCanvasSize(1200, 2300)
     h['2d_map'+str(wall)].SetWindowSize(1250, 1200)
     ut.bookCanvas(h,'cls_dis'+str(wall),'clusters distribution wall '+str(wall),cx=2,cy=5)
     h['cls_dis'+str(wall)].SetCanvasSize(1200, 2300)
     h['cls_dis'+str(wall)].SetWindowSize(1250, 1200)
     ut.bookCanvas(h,'digi_dis'+str(wall),'digit hit distribution wall '+str(wall),cx=2,cy=5)
     h['digi_dis'+str(wall)].SetCanvasSize(1200, 2300)
     h['digi_dis'+str(wall)].SetWindowSize(1250, 1200)
     ut.bookCanvas(h,'eml_map','emulsion map',1200, 1200, cx=1,cy=1)
     ut.bookHist(h, 'eml_xy', '2d emulsion map', nBins, -50, 0, nBins, 10, 60)
     for plane in range(1,6):
          ut.bookHist(h,'digi_'+str(wall)+'_'+str(plane),'2d digi wall '+str(wall)+' plane '+str(plane)+'; x [cm]; y [cm]',nBins,-50,0,nBins,10,60)
          ut.bookHist(h,'cls_'+str(wall)+'_'+str(plane),'2d cls wall '+str(wall)+' plane '+str(plane)+'; x [cm]; y [cm]',nBins,-50,0,nBins,10,60)
          ut.bookHist(h,'cls_x_'+str(wall)+'_'+str(plane),'clusters x '+str(wall)+' plane '+str(plane)+'; x [cm]',nBins,-50,0)
          ut.bookHist(h,'cls_y_'+str(wall)+'_'+str(plane),'clusters y '+str(wall)+' plane '+str(plane)+'; y [cm]',nBins,10,60)
          ut.bookHist(h,'digi_x_'+str(wall)+'_'+str(plane),'digi wall '+str(wall)+' plane '+str(plane)+'; x [cm]',nBins,-50,0)
          ut.bookHist(h,'digi_y_'+str(wall)+'_'+str(plane),'digi wall '+str(wall)+' plane '+str(plane)+'; y [cm]',nBins,10,60)
     sizeBin = 50/nBins
     #P = ROOT.TLorentzVector()
     A,B = ROOT.TVector3(),ROOT.TVector3()
     fgaus = ROOT.TF1('fgaus', 'gaus')
     for N in range(start, eventTree.GetEntries()):
          event = eventTree
          rc = event.GetEvent(N)
          motherWall = -1
          cc = 0
          hitEntries=0
          incoming_nu = event.MCTrack[0]
          outgoing_l = event.MCTrack[1]
          if incoming_nu.GetPdgCode()==14 and outgoing_l.GetPdgCode()==13:
               cc = 1
               mX = incoming_nu.GetStartX()
               mY = incoming_nu.GetStartY()
               mZ = incoming_nu.GetStartZ()
               ROOT.gGeoManager.FindNode(mX, mY, mZ)
               node = ROOT.gGeoManager.GetPath()
               #print(N, 'starting z:', incoming_nu.GetStartZ())
               #for wall in range(5):
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
                    nuE = incoming_nu.GetEnergy()
                    #t = (nuZ-mZ)/pZ
                    #nuX = mX+t*pX
                    #nuY = mY+t*pY
                    nuX = (pX/pZ)*(nuZ-mZ)+mX
                    nuY = (pY/pZ)*(nuZ-mZ)+mY
                    #print(nuX, nux, nuY, nuy)
          if motherWall == -1: continue
          if not cc: continue
          #for wall in range(5):
          #h['hit_z_wall_'+str(wall)].Reset('ICES')
          #h['cls_z_wall_'+str(wall)].Reset('ICES')
          h['eml_xy'].Reset()
          for plane in range(1, 6):
               h['cls_'+str(wall)+'_'+str(plane)].Reset('ICES')
               h['digi_'+str(wall)+'_'+str(plane)].Reset('ICES')
               h['cls_x_'+str(wall)+'_'+str(plane)].Reset('ICES')
               h['cls_y_'+str(wall)+'_'+str(plane)].Reset('ICES')
               h['digi_x_'+str(wall)+'_'+str(plane)].Reset('ICES')
               h['digi_y_'+str(wall)+'_'+str(plane)].Reset('ICES')
          if motherWall==wall:
               clusters = makeClusters(event)
               #print(N, len(clusters))
               for aCl in clusters:
                    aCl.GetPosition(A,B)
                    vertical = int(aCl.GetFirst()/100000)%10==1
                    #print(A[0], A[1], A[2], B[0], B[1], B[2])
                    #h['cls_z_wall_'+str(wall)].Fill(A[2])
                    #h['h_size'].Fill(aCl.GetN())
                    clStation = int(aCl.GetFirst()/1000000)
                    #print(N, clStation)
                    for plane in range(1, 6):
                         if plane == clStation:
                              #h['cls_'+str(wall)+'_'+str(plane)].Fill(A[0], B[1])
                              if vertical:
                                   for bin in np.arange(0, B[1]-A[1],sizeBin):
                                        h['cls_'+str(wall)+'_'+str(plane)].Fill(A[0], B[1]-bin)
                                   h['cls_x_'+str(wall)+'_'+str(plane)].Fill(A[0])
                              else:
                                   for bin in np.arange(0, B[0]-A[0], sizeBin):
                                        h['cls_'+str(wall)+'_'+str(plane)].Fill(A[0]+bin, B[1])
                                   h['cls_y_'+str(wall)+'_'+str(plane)].Fill(B[1])
               for aHit in event.Digi_ScifiHits:
                    if not aHit.isValid(): continue
                    hitEntries+=1
                    scifiDet.GetSiPMPosition(aHit.GetDetectorID(),A,B)
                    digiStation = int(aHit.GetDetectorID()/1000000)
                    #print(digiStation)
                    for plane in range(1, 6):
                         if plane == digiStation:
                              if aHit.isVertical():
                                   for bin in np.arange(0, B[1]-A[1],sizeBin):
                                        h['digi_'+str(wall)+'_'+str(plane)].Fill(A[0], B[1]-bin)
                                   h['digi_x_'+str(wall)+'_'+str(plane)].Fill(A[0])
                              else:
                                   for bin in np.arange(0, B[0]-A[0],sizeBin):
                                        h['digi_'+str(wall)+'_'+str(plane)].Fill(A[0]+bin, B[1])
                                   h['digi_y_'+str(wall)+'_'+str(plane)].Fill(B[1])
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
                    if nWall == motherWall and nPlate == 55:
                         #trackID = eHit.GetTrackID()
                         #if trackID < 0: continue
                         #if event.MCTrack[trackID].GetMotherId()==0:
                              #print('event -> ', Nevent)
                              #h['eml_z'].Fill(emlZ)
                         h['eml_xy'].Fill(emlX, emlY)
               #prevID = -1
               #for aHit in event.ScifiPoint:
               #     nHit+=1
               #     partID = aHit.GetTrackID()
               #     hitStation = aHit.station()
               #     if partID == -2: continue
               #     if partID!=prevID:
               #          hitZ = aHit.GetZ()
               #          h['hit_z_wall_'+str(wall)].Fill(hitZ)
               #          for plane in range(1, 6):
               #               if plane == hitStation:
               #                    hitX = aHit.GetX()
               #                    hitY = aHit.GetY()
               #                    #hitEloss = aHit.GetEnergyLoss()
               #                    h['hit_'+str(wall)+'_'+str(plane)].Fill(hitX, hitY)
               #                    #h['2d_sig'].Fill(hitX, hitY, hitEloss)
               #                    prevID = partID
          print('event -> ', N)
          #if nHit>0:print('number of cluster/hit', len(clusters), nHit)
          #for wall in range(5):
          #h['z_dis'].cd(1)
          #h['hit_z_wall_'+str(wall)].Draw()
          #h['z_dis'].cd(2)
          #h['cls_z_wall_'+str(wall)].Draw()
          #h['z_dis'].Update()
          for plane in range(1, 6):
               h['digi_dis'+str(wall)].cd(plane*2-1)
               if plane == wall+1:
                    maxBinX = h['digi_x_'+str(wall)+'_'+str(plane)].GetMaximumBin()
                    maxBinY = h['digi_y_'+str(wall)+'_'+str(plane)].GetMaximumBin()
                    maxValueX = h['digi_x_'+str(wall)+'_'+str(plane)].GetMaximum()
                    maxValueY = h['digi_y_'+str(wall)+'_'+str(plane)].GetMaximum()
                    maximumx=maximumy=0
                    maxBinX2 = []
                    maxBinY2 = []
                    for binx in range(1, 51):
                         valuex = h['digi_x_'+str(wall)+'_'+str(plane)].GetBinContent(binx)
                         valuey = h['digi_y_'+str(wall)+'_'+str(plane)].GetBinContent(binx)
                         if valuex == maxValueX:
                              maxBinX2.append(binx)
                         if valuey == maxValueY:
                              maxBinY2.append(binx)
                    for i in range(len(maxBinX2)):
                         maxX2 = h['digi_x_'+str(wall)+'_'+str(plane)].GetBinCenter(maxBinX2[i])
                    for i in range(len(maxBinY2)):
                         maxY2 = h['digi_y_'+str(wall)+'_'+str(plane)].GetBinCenter(maxBinY2[i])

                    maxX = h['digi_x_'+str(wall)+'_'+str(plane)].GetBinCenter(maxBinX)
                    maxY = h['digi_y_'+str(wall)+'_'+str(plane)].GetBinCenter(maxBinY)
                    print(maxX, maxY)
                    print(maxX2, maxY2)
                    #for bin in range(1, nBins+1):
                    #     #print('x', bin, h['digi_x_'+str(wall)+'_'+str(plane)].GetBinContent(bin))
                    #     if 0 < h['digi_x_'+str(wall)+'_'+str(plane)].GetBinContent(bin) < maxValueX/3:
                    #          h['digi_x_'+str(wall)+'_'+str(plane)].SetBinContent(bin, 0)
                    #     #print('y', bin, h['digi_y_'+str(wall)+'_'+str(plane)].GetBinContent(bin))
                    #     if 0 < h['digi_y_'+str(wall)+'_'+str(plane)].GetBinContent(bin) < maxValueY/3:
                    #          h['digi_y_'+str(wall)+'_'+str(plane)].SetBinContent(bin, 0)
                    print('fitx')
                    fitX = h['digi_x_'+str(wall)+'_'+str(plane)].Fit('gaus','S','',maxX-3*sizeBin,maxX+3*sizeBin)
                    fitStatusX = int(fitX)
                    print('fity')
                    fitY = h['digi_y_'+str(wall)+'_'+str(plane)].Fit('gaus','S','',maxY-3*sizeBin,maxY+3*sizeBin)
                    fitStatusY = int(fitY)
                    fitMeanX = fitX.Parameter(1)
                    fitVarX = fitX.Parameter(2)
                    fitMeanY = fitY.Parameter(1)
                    fitVarY = fitY.Parameter(2)
                    chi2X = fitX.Chi2()
                    chi2Y = fitY.Chi2()
                    ndfX = fitX.Ndf()
                    ndfY = fitY.Ndf()
                    probX = fitX.Prob()
                    probY = fitY.Prob()

                    print('ndf', ndfX, ndfY)
                    print('chi2', chi2X, chi2Y)
                    #print('chi2n', chi2X/ndfX, chi2Y/ndfY)
                    print('prob', probX, probY)
                    #fitX.Print('V')
                    #fitY.Print('V')

                    quantile = array('d',[0.5])
                    if fitStatusX!=0 or fitX.Parameter(0)>2*h['digi_x_'+str(wall)+'_'+str(plane)].GetEntries():
                         medianX = array('d',[0])
                         print('x fit not converged')
                         h['digi_x_'+str(wall)+'_'+str(plane)].GetQuantiles(1, medianX, quantile)
                         print(medianX[0])
                    if fitStatusY!=0 or fitY.Parameter(0)>2*h['digi_y_'+str(wall)+'_'+str(plane)].GetEntries():
                         medianY = array('d',[0])
                         print('y fit not converged')
                         h['digi_y_'+str(wall)+'_'+str(plane)].GetQuantiles(1, medianY, quantile)
                         print(medianY[0])
                    digiDist = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
                    digiEntries = h['digi_x_'+str(wall)+'_'+str(plane)].GetEntries()+h['digi_y_'+str(wall)+'_'+str(plane)].GetEntries()
                    print('nu vertex ({:.6f}, {:.6f})'.format(nuX, nuY))
                    print('hit - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
               h['digi_x_'+str(wall)+'_'+str(plane)].Draw()
               h['digi_dis'+str(wall)].cd(plane*2)
               h['digi_y_'+str(wall)+'_'+str(plane)].Draw()
               h['digi_dis'+str(wall)].Update() 
               h['cls_dis'+str(wall)].cd(plane*2-1)
               if plane == wall+1:
                    maxBinX = h['cls_x_'+str(wall)+'_'+str(plane)].GetMaximumBin()
                    maxBinY = h['cls_y_'+str(wall)+'_'+str(plane)].GetMaximumBin()
                    maxValueX = h['cls_x_'+str(wall)+'_'+str(plane)].GetBinContent(maxBinX)
                    maxValueY = h['cls_y_'+str(wall)+'_'+str(plane)].GetBinContent(maxBinY)
                    maxX = h['cls_x_'+str(wall)+'_'+str(plane)].GetBinCenter(maxBinX)
                    maxY = h['cls_y_'+str(wall)+'_'+str(plane)].GetBinCenter(maxBinY)
                    print(maxX, maxY)
                    #for bin in range(1, nBins+1):
                    #     #print('x', bin, h['cls_x_'+str(wall)+'_'+str(plane)].GetBinContent(bin))
                    #     if 0 < h['cls_x_'+str(wall)+'_'+str(plane)].GetBinContent(bin) < maxValueX/3:
                    #          h['cls_x_'+str(wall)+'_'+str(plane)].SetBinContent(bin, 0)
                    #     #print('y', bin, h['cls_y_'+str(wall)+'_'+str(plane)].GetBinContent(bin))
                    #     if 0 < h['cls_y_'+str(wall)+'_'+str(plane)].GetBinContent(bin) < maxValueY/3:
                    #          h['cls_y_'+str(wall)+'_'+str(plane)].SetBinContent(bin, 0)
                    fitX = h['cls_x_'+str(wall)+'_'+str(plane)].Fit('gaus','S','',maxX-3*sizeBin,maxX+3*sizeBin)
                    fitStatusX = int(fitX)
                    fitY = h['cls_y_'+str(wall)+'_'+str(plane)].Fit('gaus','S','',maxY-3*sizeBin,maxY+3*sizeBin)
                    fitStatusY = int(fitY)
                    fitMeanX = fitX.Parameter(1)
                    fitVarX = fitX.Parameter(2)
                    fitMeanY = fitY.Parameter(1)
                    fitVarY = fitY.Parameter(2)
                    chi2X = fitX.Chi2()
                    chi2Y = fitY.Chi2()
                    probX = fitX.Prob()
                    probY = fitY.Prob()
                    ndfX = fitX.Ndf()
                    ndfY = fitY.Ndf()
                    probX = fitX.Prob()
                    probY = fitY.Prob()
                    print('ndf', ndfX, ndfY)
                    print('chi2', chi2X, chi2Y)
                    #print('chi2n', chi2X/ndfX, chi2Y/ndfY)
                    print('prob', probX, probY)
                    #fitX.Print('V')
                    #fitY.Print('V')
                    quantile = array('d',[0.5])
                    if fitStatusX!=0 or fitX.Parameter(0)>2*h['cls_x_'+str(wall)+'_'+str(plane)].GetEntries():
                         medianX = array('d',[0])
                         print('x fit not converged')
                         h['cls_x_'+str(wall)+'_'+str(plane)].GetQuantiles(1, medianX, quantile)
                         fitMeanX = medianX[0]
                         print(medianX[0])
                    if fitStatusY!=0 or fitY.Parameter(0)>2*h['cls_y_'+str(wall)+'_'+str(plane)].GetEntries():
                         medianY = array('d',[0])
                         print('y fit not converged')
                         h['cls_y_'+str(wall)+'_'+str(plane)].GetQuantiles(1, medianY, quantile)
                         fitMeanY = medianY[0]
                         print(medianY[0])
                    clsDist = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
                    clsEntries = h['cls_x_'+str(wall)+'_'+str(plane)].GetEntries()+h['cls_y_'+str(wall)+'_'+str(plane)].GetEntries()
                    print('cls - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
               h['cls_x_'+str(wall)+'_'+str(plane)].Draw()
               h['cls_dis'+str(wall)].cd(plane*2)
               h['cls_y_'+str(wall)+'_'+str(plane)].Draw()
               h['2d_map'+str(wall)].cd(plane*2-1)
               h['digi_'+str(wall)+'_'+str(plane)].Draw('COLZ')
               #if plane == wall+1:
               #     maxBinXY = h['digi_'+str(wall)+'_'+str(plane)].GetMaximumBin()
               #     maxValueXY = h['digi_'+str(wall)+'_'+str(plane)].GetBinContent(maxBinXY)
               #     nx  = h['digi_'+str(wall)+'_'+str(plane)].GetXaxis().GetNbins()+2
               #     ny  = h['digi_'+str(wall)+'_'+str(plane)].GetYaxis().GetNbins()+2
               #     maxBinX2 = (maxBinXY%nx)
               #     maxBinY2 = (((maxBinXY-maxBinX2)//nx)%ny)
               #     maxX2 = h['digi_'+str(wall)+'_'+str(plane)].GetXaxis().GetBinCenter(maxBinX2)
               #     maxY2 = h['digi_'+str(wall)+'_'+str(plane)].GetYaxis().GetBinCenter(maxBinY2)
               #     print(maxX2, maxY2)
               #     #for biny in range(1, nBins+1):
               #     #     for binx in range(1, nBins+1):
               #     #          bin = h['cls_'+str(wall)+'_'+str(plane)].GetBin(binx, biny)
               #     #          #print(bin, binx, biny, h['cls_'+str(wall)+'_'+str(plane)].GetBinContent(bin))
               #     #          if h['cls_'+str(wall)+'_'+str(plane)].GetBinContent(bin) < maxValueXY/3:
               #     #               h['cls_'+str(wall)+'_'+str(plane)].SetBinContent(bin, 0)
               #     f2gaus = ROOT.TF2('xygaus', 'xygaus', maxX2-3*sizeBin,maxX2+3*sizeBin, maxY2-3*sizeBin,maxY2+3*sizeBin)
               #     f2gaus.SetParameter(1, maxX2)
               #     f2gaus.SetParameter(3, maxY2)
               #     h['digi_'+str(wall)+'_'+str(plane)].Fit('xygaus','R')
               #     fitMeanX2 = f2gaus.GetParameter(1)
               #     fitMeanY2 = f2gaus.GetParameter(3)
               #     print('2d digi - nu position = ({:.6f}, {:.6f})'.format(fitMeanX2-nuX, fitMeanY2-nuY))
               h['2d_map'+str(wall)].cd(plane*2)
               h['cls_'+str(wall)+'_'+str(plane)].Draw('COLZ')
               #if plane == wall+1:
               #     maxBinXY = h['cls_'+str(wall)+'_'+str(plane)].GetMaximumBin()
               #     maxValueXY = h['cls_'+str(wall)+'_'+str(plane)].GetBinContent(maxBinXY)
               #     nx  = h['cls_'+str(wall)+'_'+str(plane)].GetXaxis().GetNbins()+2
               #     ny  = h['cls_'+str(wall)+'_'+str(plane)].GetYaxis().GetNbins()+2
               #     maxBinX2 = (maxBinXY%nx)
               #     maxBinY2 = (((maxBinXY-maxBinX2)//nx)%ny)
               #     maxX2 = h['cls_'+str(wall)+'_'+str(plane)].GetXaxis().GetBinCenter(maxBinX2)
               #     maxY2 = h['cls_'+str(wall)+'_'+str(plane)].GetYaxis().GetBinCenter(maxBinY2)
#
               #     print(maxX2, maxY2)
               #     #for biny in range(1, nBins+1):
               #     #     for binx in range(1, nBins+1):
               #     #          bin = h['cls_'+str(wall)+'_'+str(plane)].GetBin(binx, biny)
               #     #          #print(bin, binx, biny, h['cls_'+str(wall)+'_'+str(plane)].GetBinContent(bin))
               #     #          if h['cls_'+str(wall)+'_'+str(plane)].GetBinContent(bin) < maxValueXY/3:
               #     #               h['cls_'+str(wall)+'_'+str(plane)].SetBinContent(bin, 0)
               #     f2gaus = ROOT.TF2('xygaus', 'xygaus', maxX2-3*sizeBin,maxX2+3*sizeBin, maxY2-3*sizeBin,maxY2+3*sizeBin)
               #     f2gaus.SetParameter(1, maxX2)
               #     f2gaus.SetParameter(3, maxY2)
               #     h['cls_'+str(wall)+'_'+str(plane)].Fit('xygaus','R')
               #     fitMeanX2 = f2gaus.GetParameter(1)
               #     fitMeanY2 = f2gaus.GetParameter(3)
               #     print('2d cls - nu position = ({:.6f}, {:.6f})'.format(fitMeanX2-nuX, fitMeanY2-nuY))
          h['2d_map'+str(wall)].Update()
          h['cls_dis'+str(wall)].Update()
          h['eml_map'].cd(1)
          h['eml_xy'].Draw('COLZ')
          #h['eml_map'].cd(2)
          #h['eml_z'].Draw()
          h['eml_map'].Update()
          barX = h['eml_xy'].GetMean(1)
          barY = h['eml_xy'].GetMean(2)
          print('eml x y:', barX, barY)

          print(clsDist, digiDist, clsEntries, digiEntries, hitEntries)
          #print('nu vertex ({:.6f}, {:.6f}) / cls baricenter ({:.6f}, {:.6f}) / cls 2d ({:.6f}, {:.6f})'.format(mX, mY, fitMeanX, fitMeanY, fitMeanX2, fitMeanY2))
          #print('cls 2d - nu position = ({:.6f}, {:.6f})'.format(fitMeanX2-mX, fitMeanY2-mY))
          if save:
               h['z_dis'].Print('event_'+str(N)+'_z.png')
               h['2d_map'+str(wall)].Print('event_'+str(N)+'_2d_'+str(wall)+'.png')
               h['cls_dis'+str(wall)].Print('event_'+str(N)+'_xy_'+str(wall)+'.png')
          rc = input("hit return for next event or q for quit: ")
          if rc=='q': break

def allEvents(Nev = -1):
     srX = []
     srY = []
     srD = []
     nBins = [25, 50, 100]
     nBin = 50
     sizeBin = 50/nBin
     fitRange = [3]
     #if nBin==100: fitRange = [5, 10]
     ut.bookCanvas(h,'cls_dis','clusters distribution',cx=2,cy=5)
     h['cls_dis'].SetCanvasSize(1200, 2300)
     h['cls_dis'].SetWindowSize(1250, 1200)
     for plane in range(1,6):
          ut.bookHist(h,'cls_x_'+str(plane),'clusters x plane '+str(plane)+'; x [cm]',nBin,-50,0)
          ut.bookHist(h,'cls_y_'+str(plane),'clusters y plane '+str(plane)+'; y [cm]',nBin,10,60)
     ut.bookCanvas(h,'res3','residuals3',1800, 1000, cx=3,cy=2)
     ut.bookCanvas(h,'res5','residuals5',1800, 1000, cx=3,cy=2)
     ut.bookCanvas(h,'res10','residuals10',1800, 1000, cx=3,cy=2)
     ut.bookHist(h,'res_x','x residuals; x [cm]',50,0,50)
     ut.bookHist(h,'res_vx','x variance; x [cm]',50,0,50)
     ut.bookHist(h,'res_y','y residuals; y [cm]',50,0,50)
     ut.bookHist(h,'res_vy','y variance; y [cm]',50,0,50)
     ut.bookHist(h,'res_d','distance residuals; d[cm]',50,0,50)
     ut.bookCanvas(h,'res_v_nu3','residuals dependences3',1200, 1200, cx=1,cy=2)
     ut.bookCanvas(h,'res_v_nu5','residuals dependences5',1200, 1200, cx=1,cy=2)
     ut.bookCanvas(h,'res_v_nu10','residuals dependences10',1200, 1200, cx=1,cy=2)
     h['res_v_en'] = ROOT.TGraphErrors()
     h['res_v_dp'] = ROOT.TGraphErrors()

     A,B = ROOT.TVector3(),ROOT.TVector3()
     V = ROOT.TLorentzVector()
     for binRange in fitRange:
          cccount=0
          resX = resY = resD = 0
          h['res_x'].Reset()
          h['res_y'].Reset()
          h['res_vx'].Reset()
          h['res_vy'].Reset()
          if Nev < 1: Nev = eventTree.GetEntries()
          for event in eventTree:
               if Nev<1: break
               N = eventTree.GetEntries()-Nev
               Nev=Nev-1
               motherWall = -1
               cc = 0
               for mcTrack in event.MCTrack:
                    if mcTrack.GetMotherId()==-1 and mcTrack.GetPdgCode()==14:
                         mX = mcTrack.GetStartX()
                         mY = mcTrack.GetStartY()
                         mZ = mcTrack.GetStartZ()
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
                                   pX = mcTrack.GetPx()
                                   pY = mcTrack.GetPy()
                                   pZ = mcTrack.GetPz()
                                   nuM = mcTrack.GetMass()
                                   V.SetXYZM(pX, pY, pZ, nuM)
                                   nuE = V.E()
                                   nuX = (pX/pZ)*(nuZ-mZ)+mX
                                   nuY = (pY/pZ)*(nuZ-mZ)+mY
                                   #print(nuX-mX, nuY-mY)
                    elif mcTrack.GetMotherId()==0 and mcTrack.GetPdgCode()==13:
                         cc = 1
                         break
               if motherWall == -1: continue
               if cc:
                    #print('N wall nbin fitrange', N, motherWall, nBin, binRange)
                    #print('wall', motherWall, N)
                    for plane in range(1, 6):
                         h['cls_x_'+str(plane)].Reset()
                         h['cls_y_'+str(plane)].Reset()
                    clusters = makeClusters(event)
                    #print(N, len(clusters))
                    for aCl in clusters:
                         aCl.GetPosition(A,B)
                         vertical = int(aCl.GetFirst()/100000)%10==1
                         clStation = int(aCl.GetFirst()/1000000)
                         for plane in range(1, 6):
                              if plane == clStation:
                                   if vertical: h['cls_x_'+str(plane)].Fill(A[0])
                                   else: h['cls_y_'+str(plane)].Fill(B[1])
                    #print('')
                    #print('event -> ', N)
                    for plane in range(1, 6):
                         h['cls_dis'].cd(plane*2-1)
                         if plane == motherWall+1:
                              maxBinX = h['cls_x_'+str(plane)].GetMaximumBin()
                              maxBinY = h['cls_y_'+str(plane)].GetMaximumBin()
                              maxValueX = h['cls_x_'+str(plane)].GetBinContent(maxBinX)
                              maxValueY = h['cls_y_'+str(plane)].GetBinContent(maxBinY)
                              maxX = h['cls_x_'+str(plane)].GetBinCenter(maxBinX)
                              maxY = h['cls_y_'+str(plane)].GetBinCenter(maxBinY)
                              #print('1d', maxBinX, maxBinY)
                              #for bin in range(1, nBin+1):
                              #     if h['cls_x_'+str(plane)].GetBinContent(bin) < maxValueX/3:
                              #          h['cls_x_'+str(plane)].SetBinContent(bin, 0)
                              #     if h['cls_y_'+str(plane)].GetBinContent(bin) < maxValueY/3:
                              #          h['cls_y_'+str(plane)].SetBinContent(bin, 0)
                              fitX = h['cls_x_'+str(plane)].Fit('gaus','SQ','',maxX-binRange*sizeBin,maxX+binRange*sizeBin)
                              fitStatusX = int(fitX)
                              fitY = h['cls_y_'+str(plane)].Fit('gaus','SQ','',maxY-binRange*sizeBin,maxY+binRange*sizeBin)
                              fitStatusY = int(fitY)
                              fitMeanX = fitX.Parameter(1)
                              fitVarX = fitX.Parameter(2)
                              fitMeanY = fitY.Parameter(1)
                              fitVarY = fitY.Parameter(2)
                              quantile = array('d',[0.5])
                              if fitStatusX!=0 or fitX.Parameter(0)>2*h['cls_x_'+str(plane)].GetEntries():
                                   medianX = array('d',[0])
                                   #print('x fit not converged')
                                   h['cls_x_'+str(plane)].GetQuantiles(1, medianX, quantile)
                                   fitMeanX = medianX[0]
                              if fitStatusY!=0 or fitY.Parameter(0)>2*h['cls_y_'+str(plane)].GetEntries():
                                   medianY = array('d',[0])
                                   #print('y fit not converged')
                                   h['cls_y_'+str(plane)].GetQuantiles(1, medianY, quantile)
                                   fitMeanY = medianY[0]
                              resX += abs((fitMeanX-nuX))
                              resY += abs((fitMeanY-nuY))
                              distance = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
                              resD += distance
                         h['cls_x_'+str(plane)].Draw()
                         h['cls_dis'].cd(plane*2)
                         h['cls_y_'+str(plane)].Draw()
                    h['cls_dis'].Update()
                    #stats = h['cls_x_'+str(wall+1)].FindObject('stats')
                    #stats.SetOptFit(111)
                    #print('nu vertex {:.6f} {:.6f}'.format(nuX, nuY))
                    #print('nu vertex ({:.6f}, {:.6f}) / cls baricenter ({:.6f}, {:.6f})'.format(nuX, nuY, fitMeanX, fitMeanY))
                    #print('cls - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
                    h['res_x'].Fill(abs(fitMeanX-nuX))
                    h['res_y'].Fill(abs(fitMeanY-nuY))
                    h['res_vx'].Fill(fitVarX)
                    h['res_vy'].Fill(fitVarY)
                    h['res_d'].Fill(distance)
                    #print(N, cccount, nuE, mZ, distance)
                    h['res_v_en'].SetPoint(cccount, nuE, distance)
                    h['res_v_dp'].SetPoint(cccount, mZ, distance)
                    print(binRange, N, motherWall, fitMeanX-nuX, fitMeanY-nuY, distance)
                    cccount+=1
          h['res'+str(binRange)].cd(1)
          h['res_x'].DrawNormalized()
          h['res'+str(binRange)].cd(2)
          h['res_y'].DrawNormalized()
          h['res'+str(binRange)].cd(3)
          h['res_d'].DrawNormalized()
          h['res'+str(binRange)].cd(4)
          h['res_vx'].DrawNormalized()
          h['res'+str(binRange)].cd(5)
          h['res_vy'].DrawNormalized()
          #h['res'+str(binRange)].Update()
          h['res_v_nu'+str(binRange)].cd(1)
          h['res_v_en'].SetTitle("neutrino energy")
          h['res_v_en'].GetXaxis().SetTitle("energy [MeV]")
          h['res_v_en'].GetYaxis().SetTitle("residual [cm]")
          h['res_v_en'].SetMarkerStyle(20)
          h['res_v_en'].SetMarkerSize(0.5)
          h['res_v_en'].Draw('AP')
          h['res_v_nu'+str(binRange)].cd(2)
          h['res_v_dp'].SetTitle("neutrino interaction depth")
          h['res_v_dp'].GetXaxis().SetTitle("z [cm]")
          h['res_v_dp'].GetYaxis().SetTitle("residual [cm]")
          h['res_v_dp'].SetMarkerStyle(20)
          h['res_v_dp'].SetMarkerSize(0.5)
          h['res_v_dp'].Draw('AP')
          #h['res_v_nu'].Update()
          #h['res'+str(binRange)].Print('/home/fabio/Immagini/220725/residuals_'+str(nBin)+'_'+str(binRange)+'.png')
          #h['res_v_nu'+str(binRange)].Print('/home/fabio/Immagini/220725/nu_'+str(nBin)+'_'+str(binRange)+'.png')

          srX.append(resX)
          srY.append(resY)
          srD.append(resD)
     print('srX', srX)
     print('srY', srY)
     print('srD', srD)

def allDigit(Nev = -1):
     srX = []
     srY = []
     srD = []
     nBins = [25, 50, 100]
     nBin = 50
     sizeBin = 50/nBin
     fitRange = [3]
     #if nBin==100: fitRange = [5, 10]
     ut.bookCanvas(h,'digi_dis','digi_htis distribution',cx=2,cy=5)
     h['digi_dis'].SetCanvasSize(1200, 2300)
     h['digi_dis'].SetWindowSize(1250, 1200)
     for plane in range(1,6):
          ut.bookHist(h,'digi_x_'+str(plane),'digi_hits x plane '+str(plane)+'; x [cm]',nBin,-50,0)
          ut.bookHist(h,'digi_y_'+str(plane),'digi_hits y plane '+str(plane)+'; y [cm]',nBin,10,60)
     ut.bookCanvas(h,'res3','residuals3',1800, 1000, cx=3,cy=2)
     #ut.bookCanvas(h,'res5','residuals5',1800, 1000, cx=3,cy=2)
     #ut.bookCanvas(h,'res10','residuals10',1800, 1000, cx=3,cy=2)
     ut.bookHist(h,'res_x','x residuals; x [cm]',50,0,50)
     ut.bookHist(h,'res_vx','x variance; x [cm]',50,0,50)
     ut.bookHist(h,'res_y','y residuals; y [cm]',50,0,50)
     ut.bookHist(h,'res_vy','y variance; y [cm]',50,0,50)
     ut.bookHist(h,'res_d','distance residuals; [cm]',50,0,50)
     ut.bookCanvas(h,'res_v_nu','residuals dependences',1200, 1200, cx=1,cy=2)
     h['res_v_en'] = ROOT.TGraphErrors()
     h['res_v_dp'] = ROOT.TGraphErrors()
     A,B = ROOT.TVector3(),ROOT.TVector3()
     V = ROOT.TLorentzVector()
     for binRange in fitRange:
          cccount=0
          resX = resY = resD = 0
          h['res_x'].Reset()
          h['res_y'].Reset()
          h['res_vx'].Reset()
          h['res_vy'].Reset()
          if Nev < 1: Nev = eventTree.GetEntries()
          for event in eventTree:
               if Nev<1: break
               N = eventTree.GetEntries()-Nev
               Nev=Nev-1
               motherWall = -1
               cc = 0
               for mcTrack in event.MCTrack:
                    if mcTrack.GetMotherId()==-1 and mcTrack.GetPdgCode()==14:
                         mX = mcTrack.GetStartX()
                         mY = mcTrack.GetStartY()
                         mZ = mcTrack.GetStartZ()
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
                                   pX = mcTrack.GetPx()
                                   pY = mcTrack.GetPy()
                                   pZ = mcTrack.GetPz()
                                   nuM = mcTrack.GetMass()
                                   V.SetXYZM(pX, pY, pZ, nuM)
                                   nuE = V.E()
                                   nuX = (pX/pZ)*(nuZ-mZ)+mX
                                   nuY = (pY/pZ)*(nuZ-mZ)+mY
                    elif mcTrack.GetMotherId()==0 and mcTrack.GetPdgCode()==13:
                         cc = 1
                         break
               if motherWall == -1: continue
               if cc:
                    #print('N wall nbin fitrange', N, motherWall, nBin, binRange)
                    #print('wall', motherWall, N)
                    for plane in range(1, 6):
                         h['digi_x_'+str(plane)].Reset()
                         h['digi_y_'+str(plane)].Reset()
                    for aHit in event.Digi_ScifiHits:
                         if not aHit.isValid(): continue
                         scifiDet.GetSiPMPosition(aHit.GetDetectorID(),A,B)
                         digiStation = int(aHit.GetDetectorID()/1000000)
                         for plane in range(1, 6):
                              if plane == digiStation:
                                   if aHit.isVertical():
                                        #for bin in np.arange(0, B[1]-A[1],sizeBin):
                                        #     h['digi_'+str(plane)].Fill(A[0], B[1]-bin)
                                        h['digi_x_'+str(plane)].Fill(A[0])
                                   else:
                                        #for bin in np.arange(0, B[0]-A[0],sizeBin):
                                        #     h['digi_'+str(plane)].Fill(A[0]+bin, B[1])
                                        h['digi_y_'+str(plane)].Fill(B[1])
                    #print('')
                    #print('event -> ', N)
                    for plane in range(1, 6):
                         h['digi_dis'].cd(plane*2-1)
                         if plane == motherWall+1:
                              maxBinX = h['digi_x_'+str(plane)].GetMaximumBin()
                              maxBinY = h['digi_y_'+str(plane)].GetMaximumBin()
                              maxValueX = h['digi_x_'+str(plane)].GetBinContent(maxBinX)
                              maxValueY = h['digi_y_'+str(plane)].GetBinContent(maxBinY)
                              maxX = h['digi_x_'+str(plane)].GetBinCenter(maxBinX)
                              maxY = h['digi_y_'+str(plane)].GetBinCenter(maxBinY)
                              #print('1d', maxBinX, maxBinY)
                              #for bin in range(1, nBin+1):
                              #     if h['digi_x_'+str(plane)].GetBinContent(bin) < maxValueX/3:
                              #          h['digi_x_'+str(plane)].SetBinContent(bin, 0)
                              #     if h['digi_y_'+str(plane)].GetBinContent(bin) < maxValueY/3:
                              #          h['digi_y_'+str(plane)].SetBinContent(bin, 0)
                              fitX = h['digi_x_'+str(plane)].Fit('gaus','SQ','',maxX-binRange*sizeBin,maxX+binRange*sizeBin)
                              fitStatusX = int(fitX)
                              fitY = h['digi_y_'+str(plane)].Fit('gaus','SQ','',maxY-binRange*sizeBin,maxY+binRange*sizeBin)
                              fitStatusY = int(fitY)
                              fitMeanX = fitX.Parameter(1)
                              fitVarX = fitX.Parameter(2)
                              fitMeanY = fitY.Parameter(1)
                              fitVarY = fitY.Parameter(2)
                              quantile = array('d',[0.5])
                              if fitStatusX!=0 or fitX.Parameter(0)>2*h['digi_x_'+str(plane)].GetEntries():
                                   medianX = array('d',[0])
                                   #print('x fit not converged')
                                   h['digi_x_'+str(plane)].GetQuantiles(1, medianX, quantile)
                                   fitMeanX = medianX[0]
                              if fitStatusY!=0 or fitY.Parameter(0)>2*h['digi_y_'+str(plane)].GetEntries():
                                   medianY = array('d',[0])
                                   #print('y fit not converged')
                                   h['digi_y_'+str(plane)].GetQuantiles(1, medianY, quantile)
                                   fitMeanY = medianY[0]
                              resX += abs((fitMeanX-nuX))
                              resY += abs((fitMeanY-nuY))
                              distance = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
                              resD += distance
                         h['digi_x_'+str(plane)].Draw()
                         h['digi_dis'].cd(plane*2)
                         h['digi_y_'+str(plane)].Draw()
                    h['digi_dis'].Update()
                    #stats = h['digi_x_'+str(wall+1)].FindObject('stats')
                    #stats.SetOptFit(111)
                    #print('nu vertex {:.6f} {:.6f}'.format(nuX, nuY))
                    #print('nu vertex ({:.6f}, {:.6f}) / digi baricenter ({:.6f}, {:.6f})'.format(nuX, nuY, fitMeanX, fitMeanY))
                    #print('digi - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
                    h['res_x'].Fill(abs(fitMeanX-nuX))
                    h['res_y'].Fill(abs(fitMeanY-nuY))
                    h['res_vx'].Fill(fitVarX)
                    h['res_vy'].Fill(fitVarY)
                    h['res_d'].Fill(distance)
                    #print(N, cccount, nuE, mZ, distance)
                    h['res_v_en'].SetPoint(cccount, nuE, distance)
                    h['res_v_dp'].SetPoint(cccount, mZ, distance)
                    print(N, motherWall, fitMeanX-nuX, fitMeanY-nuY, distance)
                    cccount+=1
          h['res'+str(binRange)].cd(1)
          h['res_x'].DrawNormalized()
          h['res'+str(binRange)].cd(2)
          h['res_y'].DrawNormalized()
          h['res'+str(binRange)].cd(3)
          h['res_d'].DrawNormalized()
          h['res'+str(binRange)].cd(4)
          h['res_vx'].DrawNormalized()
          h['res'+str(binRange)].cd(5)
          h['res_vy'].DrawNormalized()
          #h['res'+str(binRange)].Update()
          h['res_v_nu'].cd(1)
          h['res_v_en'].SetTitle("neutrino energy")
          h['res_v_en'].GetXaxis().SetTitle("energy [MeV]")
          h['res_v_en'].GetYaxis().SetTitle("residual [cm]")
          h['res_v_en'].SetMarkerStyle(20)
          h['res_v_en'].SetMarkerSize(0.5)
          h['res_v_en'].Draw('AP')
          h['res_v_nu'].cd(2)
          h['res_v_dp'].SetTitle("neutrino interaction depth")
          h['res_v_dp'].GetXaxis().SetTitle("z [cm]")
          h['res_v_dp'].GetYaxis().SetTitle("residual [cm]")
          h['res_v_dp'].SetMarkerStyle(20)
          h['res_v_dp'].SetMarkerSize(0.5)
          h['res_v_dp'].Draw('AP')
          #h['res_v_nu'].Update()
          #h['res'+str(binRange)].Print('/home/fabio/Immagini/220721/residuals_'+str(nBin)+'_'+str(binRange)+'_ext_thr.png')
          srX.append(resX)
          srY.append(resY)
          srD.append(resD)
     print('srX', srX)
     print('srY', srY)
     print('srD', srD)

def allall(Nev = -1):
     nBin = 50
     sizeBin = 50/nBin
     fitRange = [3]
     
     ut.bookCanvas(h,'cls_dis','clusters distribution',cx=2,cy=5)
     h['cls_dis'].SetCanvasSize(1200, 2300)
     h['cls_dis'].SetWindowSize(1250, 1200)
     for plane in range(1,6):
          ut.bookHist(h,'cls_x_'+str(plane),'clusters x plane '+str(plane)+'; x [cm]',nBin,-50,0)
          ut.bookHist(h,'cls_y_'+str(plane),'clusters y plane '+str(plane)+'; y [cm]',nBin,10,60)
     ut.bookCanvas(h,'digi_dis','digi_hits distribution',cx=2,cy=5)
     h['digi_dis'].SetCanvasSize(1200, 2300)
     h['digi_dis'].SetWindowSize(1250, 1200)
     for plane in range(1,6):
          ut.bookHist(h,'digi_x_'+str(plane),'digi_hits x plane '+str(plane)+'; x [cm]',nBin,-50,0)
          ut.bookHist(h,'digi_y_'+str(plane),'digi_hits y plane '+str(plane)+'; y [cm]',nBin,10,60)
     ut.bookCanvas(h,'cls_v_digi','better fit',1200, 1200, cx=2,cy=2)
     ut.bookHist(h,'digi_fitg','digi_hits good fit; #hits',45,0,900)
     ut.bookHist(h,'cls_fitg','cls good fit;#clusters',60,0,120)
     ut.bookHist(h,'digi_fitb','digi_hits bad fit;#hits',45,0,900)
     ut.bookHist(h,'cls_fitb','cls bad fit;#clusters',60,0,120)
     ut.bookHist(h,'clh_fitg','hits in cls good fit;#hits',20,0,20)
     ut.bookHist(h,'clh_fitb','hits in cls bad fit;#hits',20,0,20)
     ut.bookHist(h,'diffc','digi-cls (cls); d[cm]',50,-25,25)
     ut.bookHist(h,'diffd','digi-cls (digi); d[cm]',50,-25,25)

     A,B = ROOT.TVector3(),ROOT.TVector3()
     V = ROOT.TLorentzVector()
     for binRange in fitRange:
          clse = []
          digie = []
          cccount=0
          fitcls=0
          fitdigis=0
          if Nev < 1: Nev = eventTree.GetEntries()
          for event in eventTree:
               if Nev<1: break
               N = eventTree.GetEntries()-Nev
               Nev=Nev-1
               motherWall = -1
               cc = 0
               digicount=0
               incoming_nu = event.MCTrack[0]
               outgoing_l = event.MCTrack[1]
               clsHits = []
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
                              nuE = incoming_nu.GetEnergy()
                              nuX = (pX/pZ)*(nuZ-mZ)+mX
                              nuY = (pY/pZ)*(nuZ-mZ)+mY
               if motherWall == -1: continue
               if not cc: continue
               #print('N wall nbin fitrange', N, motherWall, nBin, binRange)
               #print('wall', motherWall, N)
               for plane in range(1, 6):
                    h['cls_x_'+str(plane)].Reset()
                    h['cls_y_'+str(plane)].Reset()
                    h['digi_x_'+str(plane)].Reset()
                    h['digi_y_'+str(plane)].Reset()
               clusters = makeClusters(event)
               #print(N, len(clusters))
               for aHit in event.Digi_ScifiHits:
                    if not aHit.isValid(): continue
                    digicount+=1
                    scifiDet.GetSiPMPosition(aHit.GetDetectorID(),A,B)
                    digiStation = int(aHit.GetDetectorID()/1000000)
                    for plane in range(1, 6):
                         if plane == digiStation:
                              if aHit.isVertical(): h['digi_x_'+str(plane)].Fill(A[0])
                              else: h['digi_y_'+str(plane)].Fill(B[1])
               for aCl in clusters:
                    aCl.GetPosition(A,B)
                    vertical = int(aCl.GetFirst()/100000)%10==1
                    clStation = int(aCl.GetFirst()/1000000)
                    clsHits.append(aCl.GetN())
                    for plane in range(1, 6):
                         if plane == clStation:
                              if vertical: h['cls_x_'+str(plane)].Fill(A[0])
                              else: h['cls_y_'+str(plane)].Fill(B[1])
               print('event -> ', N)
               if h['digi_y_'+str(motherWall+1)].GetEntries()==0 or h['digi_x_'+str(motherWall+1)].GetEntries()==0:
                    continue
               if h['cls_y_'+str(motherWall+1)].GetEntries()==0 or h['cls_x_'+str(motherWall+1)].GetEntries()==0:
                    continue
               if digicount > 100: continue
               for plane in range(1, 6):
                    h['cls_dis'].cd(plane*2-1)
                    if plane == motherWall+1:
                         maxBinX = h['cls_x_'+str(plane)].GetMaximumBin()
                         maxBinY = h['cls_y_'+str(plane)].GetMaximumBin()
                         maxValueX = h['cls_x_'+str(plane)].GetBinContent(maxBinX)
                         maxValueY = h['cls_y_'+str(plane)].GetBinContent(maxBinY)
                         maxX = h['cls_x_'+str(plane)].GetBinCenter(maxBinX)
                         maxY = h['cls_y_'+str(plane)].GetBinCenter(maxBinY)
                         fitX = h['cls_x_'+str(plane)].Fit('gaus','SQ','',maxX-binRange*sizeBin,maxX+binRange*sizeBin)
                         fitStatusX = int(fitX)
                         fitY = h['cls_y_'+str(plane)].Fit('gaus','SQ','',maxY-binRange*sizeBin,maxY+binRange*sizeBin)
                         fitStatusY = int(fitY)
                         fitMeanX = fitX.Parameter(1)
                         fitVarX = fitX.Parameter(2)
                         fitMeanY = fitY.Parameter(1)
                         fitVarY = fitY.Parameter(2)
                         quantile = array('d',[0.5])
                         if fitStatusX!=0 or fitX.Parameter(0)>2*h['cls_x_'+str(plane)].GetEntries():
                              medianX = array('d',[0])
                              #print('x fit not converged')
                              h['cls_x_'+str(plane)].GetQuantiles(1, medianX, quantile)
                              fitMeanX = medianX[0]
                         if fitStatusY!=0 or fitY.Parameter(0)>2*h['cls_y_'+str(plane)].GetEntries():
                              medianY = array('d',[0])
                              #print('y fit not converged')
                              h['cls_y_'+str(plane)].GetQuantiles(1, medianY, quantile)
                              fitMeanY = medianY[0]
                         clsDist = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
                         clsEntries = h['cls_x_'+str(plane)].GetEntries()+h['cls_y_'+str(plane)].GetEntries()
                         clse.append(clsEntries)
                         print('nu vertex ({:.6f}, {:.6f})'.format(nuX, nuY))
                         print('cls - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
                    h['cls_x_'+str(plane)].Draw()
                    h['cls_dis'].cd(plane*2)
                    h['cls_y_'+str(plane)].Draw()
                    h['digi_dis'].cd(plane*2-1)
                    if plane == motherWall+1:
                         maxBinX = h['digi_x_'+str(plane)].GetMaximumBin()
                         maxBinY = h['digi_y_'+str(plane)].GetMaximumBin()
                         maxValueX = h['digi_x_'+str(plane)].GetBinContent(maxBinX)
                         maxValueY = h['digi_y_'+str(plane)].GetBinContent(maxBinY)
                         maxX = h['digi_x_'+str(plane)].GetBinCenter(maxBinX)
                         maxY = h['digi_y_'+str(plane)].GetBinCenter(maxBinY)
                         fitX = h['digi_x_'+str(plane)].Fit('gaus','SQ','',maxX-3*sizeBin,maxX+3*sizeBin)
                         fitStatusX = int(fitX)
                         fitY = h['digi_y_'+str(plane)].Fit('gaus','SQ','',maxY-3*sizeBin,maxY+3*sizeBin)
                         fitStatusY = int(fitY)
                         fitMeanX = fitX.Parameter(1)
                         fitVarX = fitX.Parameter(2)
                         fitMeanY = fitY.Parameter(1)
                         fitVarY = fitY.Parameter(2)
                         quantile = array('d',[0.5])
                         if fitStatusX!=0 or fitX.Parameter(0)>2*h['digi_x_'+str(plane)].GetEntries():
                              medianX = array('d',[0])
                              h['digi_x_'+str(plane)].GetQuantiles(1, medianX, quantile)
                         if fitStatusY!=0 or fitY.Parameter(0)>2*h['digi_y_'+str(plane)].GetEntries():
                              medianY = array('d',[0])
                              h['digi_y_'+str(plane)].GetQuantiles(1, medianY, quantile)
                         digiDist = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
                         digiEntries = h['digi_x_'+str(plane)].GetEntries()+h['digi_y_'+str(plane)].GetEntries()
                         digie.append(digiEntries)
                         print('hit - nu position = ({:.6f}, {:.6f})'.format(fitMeanX-nuX, fitMeanY-nuY))
                    h['digi_x_'+str(plane)].Draw()
                    h['digi_dis'].cd(plane*2)
                    h['digi_y_'+str(plane)].Draw()
               if digiEntries<80 and clsEntries<30: pass
               h['cls_dis'].Update()
               h['digi_dis'].Update()
               print(clsDist, digiDist, clsEntries, digiEntries)
               if clsDist<digiDist:
                    fitcls+=1
                    h['diffc'].Fill(digiDist-clsDist)
                    h['cls_fitg'].Fill(clsEntries)
                    h['digi_fitb'].Fill(digiEntries)
                    for nhit in clsHits:
                         h['clh_fitg'].Fill(nhit)
                    #print('cls->', clsEntries, digiEntries)
               else: 
                    fitdigis+=1
                    h['diffd'].Fill(digiDist-clsDist)
                    h['cls_fitb'].Fill(clsEntries)
                    h['digi_fitg'].Fill(digiEntries)
                    for nhit in clsHits:
                         h['clh_fitb'].Fill(nhit)
                    #print('digi->', digiEntries, clsEntries)
               cccount+=1
          #print('cls/digi', max(clse), max(digie))
          #print(max(clsHits))
          print('cls/digis', cccount, fitcls/cccount, fitdigis/cccount)
          h['cls_fitg'].SetLineColor(1)
          h['cls_fitb'].SetLineColor(2)
          h['digi_fitg'].SetLineColor(1)
          h['digi_fitb'].SetLineColor(2)
          h['clh_fitg'].SetLineColor(1)
          h['clh_fitb'].SetLineColor(2)
          h['diffc'].SetLineColor(1)
          h['diffd'].SetLineColor(2)
          h['cls_fitg'].SetStats(False)
          h['cls_fitb'].SetStats(False)
          h['digi_fitg'].SetStats(False)
          h['digi_fitb'].SetStats(False)
          h['clh_fitg'].SetStats(False)
          h['clh_fitb'].SetStats(False)
          h['diffc'].SetStats(False)
          h['diffd'].SetStats(False)
          h['cls_v_digi'].cd(1)
          h['cls_fitb'].Draw()
          h['cls_fitg'].Draw('SAME')
          h['cls_v_digi'].cd(2)
          h['digi_fitb'].Draw()
          h['digi_fitg'].Draw('SAME')
          h['cls_v_digi'].cd(3)
          h['clh_fitb'].Draw()
          h['clh_fitg'].Draw('SAME')
          h['cls_v_digi'].cd(4)
          h['diffc'].Draw()
          h['diffd'].Draw('SAME')
          h['cls_v_digi'].cd(1).BuildLegend(0.6, 0.7, 0.9, 0.9)
          h['cls_v_digi'].cd(2).BuildLegend(0.6, 0.7, 0.9, 0.9)
          h['cls_v_digi'].cd(3).BuildLegend(0.6, 0.7, 0.9, 0.9)
          h['cls_v_digi'].cd(4).BuildLegend(0.6, 0.7, 0.9, 0.9)

def sigPos(Nev = -1):
     nBin = 50
     sizeBin = 50/nBin
     binRange = 3
     ut.bookCanvas(h,'pos','plane positions',1200, 500, cx=2,cy=1)
     ut.bookHist(h,'posX','x position ; x [cm]',nBin,-50,0)
     ut.bookHist(h,'posY','y position ; y [cm]',nBin,10,60)
     ut.bookCanvas(h,'res','residuals',1800, 1000, cx=3,cy=2)
     ut.bookHist(h,'resX','x residuals; x [cm]',50,0,5)
     ut.bookHist(h,'resVx','x variance; x [cm]',50,0,5)
     ut.bookHist(h,'resY','y residuals; y [cm]',50,0,5)
     ut.bookHist(h,'resVy','y variance; y [cm]',50,0,5)
     ut.bookHist(h,'resD','distance residuals; d[cm]',50,0,5)
     ut.bookCanvas(h,'res_v_fit','residuals chi2',1200, 1200, cx=2,cy=2)
     ut.bookHist(h,'reschi2_x',' ;resx;chi2/ndf',50,0,20, 50, 0, 20)
     ut.bookHist(h,'resp_x',' ;resx;p',50,0,20, 50, 0, 1)
     ut.bookHist(h,'reschi2_y',' ;resy;chi2/ndf',50,0,20, 50, 0, 20)
     ut.bookHist(h,'resp_y',' ;resy;p',50,0,20, 40, 0, 1)
     ut.bookCanvas(h,'zerodf', 'zerodf', 800, 800, 1, 1)
     ut.bookHist(h, 'no_df', 'zero df', 50, 0, 25)
     ndf = 0
     A, B = ROOT.TVector3(), ROOT.TVector3()
     ccCount=0
     srX = srY = srD = 0
     fitcls=0
     fitdigis=0
     if Nev < 0: Nev = eventTree.GetEntries()
     for event in eventTree:
          if Nev<1: break
          Nevent = eventTree.GetEntries()-Nev
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
                         nuE = incoming_nu.GetEnergy()
                         nuX = (pX/pZ)*(nuZ-mZ)+mX
                         nuY = (pY/pZ)*(nuZ-mZ)+mY
          if motherWall == -1: continue
          if not cc: continue
          #print('N wall nbin fitrange', N, motherWall, nBin, binRange)
          #print('wall', motherWall, N)
          h['posX'].Reset()
          h['posY'].Reset()
          #digiSciFi = ROOT.TClonesArray("sndScifiHit")
          digiSciFi = event.Digi_ScifiHits
          hitDict = {}
          for k in range(digiSciFi.GetEntries()):
               d = digiSciFi[k]
               if not d.isValid(): continue 
               hitStation = int(d.GetDetectorID()/1000000)
               if hitStation == motherWall+1: hitDict[d.GetDetectorID()] = k
          hitList = list(hitDict.keys())
          if len(hitList) > 80:
               whichFit = 'digi'
               fitdigis+=1
               for aHit in hitList:
                    scifiDet.GetSiPMPosition(aHit,A,B)
                    #digiStation = int(aHit/1000000)
                    #if digiStation == motherWall+1:
                    vertical = int(aHit/100000)%10==1
                    if vertical: h['posX'].Fill(A[0])
                    else: h['posY'].Fill(B[1])
          elif len(hitList) > 0:
               whichFit = 'cls'
               fitcls+=1
               clusScifi = ROOT.TClonesArray("sndCluster")
               index = 0
               hitList.sort()
               tmp = [ hitList[0] ]
               cprev = hitList[0]
               ncl = 0
               last = len(hitList)-1
               hitlist = ROOT.std.vector("sndScifiHit*")()
               for i in range(len(hitList)):
                    if i==0 and len(hitList)>1: continue
                    c=hitList[i]
                    neighbour = False
                    if (c-cprev)==1:    # does not account for neighbours across sipms
                         neighbour = True
                         tmp.append(c)
                    if not neighbour  or c==hitList[last]:
                         first = tmp[0]
                         N = len(tmp)
                         hitlist.clear()
                         for aHit in tmp: hitlist.push_back( digiSciFi[hitDict[aHit]])
                         aCluster = ROOT.sndCluster(first,N,hitlist,scifiDet,False)
                         nhit = aCluster.GetN()
                         if  clusScifi.GetSize() == index: clusScifi.Expand(index+10)
                         clusScifi[index]=aCluster
                         index+=1
                         if c!=hitList[last]:
                              ncl+=1
                              tmp = [c]
                         elif not neighbour :   # save last channel
                              hitlist.clear()
                              hitlist.push_back( digiSciFi[hitDict[c]])
                              aCluster = ROOT.sndCluster(c,1,hitlist,scifiDet,False)
                              clusScifi[index]=aCluster
                              index+=1
                    cprev = c
               for aCl in clusScifi:
                    aCl.GetPosition(A,B)
                    #clStation = int(aCl.GetFirst()/1000000)
                    #if clStation == motherWall+1:
                    vertical = int(aCl.GetFirst()/100000)%10==1
                    if vertical: h['posX'].Fill(A[0])
                    else: h['posY'].Fill(B[1])
          #if h['posX'].GetEntries()==0 or h['posY'].GetEntries()==0:
          #     continue
          #print('event -> ', Nevent)
          h['pos'].cd(1)
          maxBinX = h['posX'].GetMaximumBin()
          maxBinY = h['posY'].GetMaximumBin()
          #maxValueX = h['posX'].GetBinContent(maxBinX)
          #maxValueY = h['posY'].GetBinContent(maxBinY)
          maxX = h['posX'].GetBinCenter(maxBinX)
          maxY = h['posY'].GetBinCenter(maxBinY)
          fitX = h['posX'].Fit('gaus','SQ','',maxX-binRange*sizeBin,maxX+binRange*sizeBin)
          fitStatusX = int(fitX)
          fitY = h['posY'].Fit('gaus','SQ','',maxY-binRange*sizeBin,maxY+binRange*sizeBin)
          fitStatusY = int(fitY)
          fitMeanX = fitX.Parameter(1)
          fitVarX = fitX.Parameter(2)
          fitMeanY = fitY.Parameter(1)
          fitVarY = fitY.Parameter(2)
          chi2X = fitX.Chi2()
          chi2Y = fitY.Chi2()
          probX = fitX.Prob()
          probY = fitY.Prob()
          ndfX = fitX.Ndf()
          ndfY = fitY.Ndf()
          quantile = array('d',[0.5])
          if fitStatusX!=0 or fitX.Parameter(0)>2*h['posX'].GetEntries():
               medianX = array('d',[0])
               h['posX'].GetQuantiles(1, medianX, quantile)
               fitMeanX = medianX[0]
          if fitStatusY!=0 or fitY.Parameter(0)>2*h['posY'].GetEntries():
               medianY = array('d',[0])
               h['posY'].GetQuantiles(1, medianY, quantile)
               fitMeanY = medianY[0]
          dist = ROOT.TMath.Sqrt((fitMeanX-nuX)**2+(fitMeanY-nuY)**2)
          srX += abs((fitMeanX-nuX))
          srY += abs((fitMeanY-nuY))
          srD += dist
          print('event {} wall {} res ({:.6f}, {:.6f} -> {:.6f})  '.format(Nevent,motherWall, fitMeanX-nuX, fitMeanY-nuY, dist)+whichFit)
          h['posX'].Draw()
          h['pos'].cd(2)
          h['posY'].Draw()
          h['pos'].Update()
          h['resX'].Fill(abs(fitMeanX-nuX))
          h['resY'].Fill(abs(fitMeanY-nuY))
          h['resVx'].Fill(fitVarX)
          h['resVy'].Fill(fitVarY)
          h['resD'].Fill(dist)
          if ndfX>0:
               h['reschi2_x'].Fill(abs(fitMeanX-nuX), chi2X/ndfX)
          else:
               ndf+=1
               h['no_df'].Fill(abs(fitMeanX-nuX))
          if ndfY > 0 :
               h['reschi2_y'].Fill(abs(fitMeanY-nuY), chi2Y/ndfY)
          else:
               ndf+=1
               h['no_df'].Fill(abs(fitMeanY-nuY))
          h['resp_x'].Fill(abs(fitMeanX-nuX), probX)
          h['resp_y'].Fill(abs(fitMeanY-nuY), probY)
          ccCount+=1
     for i, iRes in enumerate(['X', 'Y', 'D', 'Vx', 'Vy']):
          h['res'].cd(i+1)
          h['res'+str(iRes)].DrawNormalized()
     print('cc events', ccCount, fitdigis/ccCount, fitcls/ccCount)
     print('srX', srX)
     print('srY', srY)
     print('srD', srD)
     print('ndf', ndf/2/ccCount)
     h['res_v_fit'].cd(1)
     h['reschi2_x'].Draw('COLZ')
     h['res_v_fit'].cd(2)
     h['resp_x'].Draw('COLZ')
     h['res_v_fit'].cd(3)
     h['reschi2_y'].Draw('COLZ')
     h['res_v_fit'].cd(4)
     h['resp_y'].Draw('COLZ')
     h['res_v_fit'].Update()
     h['zerodf'].cd()
     h['no_df'].Draw()

def charts():
     nx = 5+5+5
     ut.bookCanvas(h,'charts','',1000,500,1,1) 
     ut.bookHist(h,'bar_chart','',nx,0,nx)
     ut.bookHist(h,'bar_chart_mu','',nx,0,nx)
     labels = ['Wall_0', 'Wall_1', 'Wall_2', 'Wall_3', 'Wall_4',
               'Border_0', 'Border_1', 'Border_2', 'Border_3', 'Border_4',
               'Scifi_1', 'Scifi_2', 'Scifi_3', 'Scifi_4', 'Scifi_5'] 
     for event in eventTree:
          for mcTrack in event.MCTrack:
               if mcTrack.GetMotherId()==-1:
                    mX = mcTrack.GetStartX()
                    mY = mcTrack.GetStartY()
                    mZ = mcTrack.GetStartZ()
                    ROOT.gGeoManager.FindNode(mX, mY, mZ)
                    node = ROOT.gGeoManager.GetPath()
                    for i in range(5):
                         if 'Wall_'+str(i) in node:
                              h['bar_chart'].Fill(i)
                              if mcTrack.GetPdgCode()==14:
                                   h['bar_chart_mu'].Fill(i)
                         elif 'Wallborder_'+str(i) in node:
                              h['bar_chart'].Fill(i+5)
                              if mcTrack.GetPdgCode()==14:
                                   h['bar_chart_mu'].Fill(i+5)
                         elif 'ScifiVolume'+str(i+1) in node:
                              h['bar_chart'].Fill(i+10)
                              if mcTrack.GetPdgCode()==14:
                                   h['bar_chart_mu'].Fill(i+10)
     for i in range(nx):
          h['bar_chart'].GetXaxis().SetBinLabel(i+1, labels[i])
     h['bar_chart'].SetFillColor(4)
     h['bar_chart'].SetBarWidth(0.4)
     h['bar_chart'].SetBarOffset(0.1)
     #h['bar_chart'].SetStats(0)
     
     h['bar_chart_mu'].SetFillColor(38)
     h['bar_chart_mu'].SetBarWidth(0.4)
     h['bar_chart_mu'].SetBarOffset(0.5)
     #h['bar_chart_mu'].SetStats(0)
     
     h['bar_chart'].Draw("b")
     h['bar_chart_mu'].Draw("b&&same")

     legend = ROOT.TLegend(0.6, 0.55, 0.8, 0.7)
     legend.AddEntry(h['bar_chart'],"All particles","F")
     legend.AddEntry(h['bar_chart_mu'],"Muonic neutrinos","F")
     legend.Draw()

def makeClusters(event):
     scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
     clusScifi = ROOT.TClonesArray("sndCluster")
     #digiSciFi = ROOT.TClonesArray("sndScifiHit")
     digiSciFi = event.Digi_ScifiHits
     index = 0
     hitDict = {}
     for k in range(digiSciFi.GetEntries()):
          d = digiSciFi[k]
          if not d.isValid(): continue 
          hitDict[d.GetDetectorID()] = k
     hitList = list(hitDict.keys())
     if len(hitList)>0:
          hitList.sort()
          tmp = [ hitList[0] ]
          cprev = hitList[0]
          ncl = 0
          last = len(hitList)-1
          hitlist = ROOT.std.vector("sndScifiHit*")()
          for i in range(len(hitList)):
               if i==0 and len(hitList)>1: continue
               c=hitList[i]
               neighbour = False
               if (c-cprev)==1:    # does not account for neighbours across sipms
                    neighbour = True
                    tmp.append(c)
               if not neighbour  or c==hitList[last]:
                    first = tmp[0]
                    N = len(tmp)
                    hitlist.clear()
                    for aHit in tmp: hitlist.push_back( digiSciFi[hitDict[aHit]])
                    aCluster = ROOT.sndCluster(first,N,hitlist,scifiDet,False)
                    nhit = aCluster.GetN()
                    if  clusScifi.GetSize() == index: clusScifi.Expand(index+10)
                    clusScifi[index]=aCluster
                    index+=1
                    if c!=hitList[last]:
                         ncl+=1
                         tmp = [c]
                    elif not neighbour :   # save last channel
                         hitlist.clear()
                         hitlist.push_back( digiSciFi[hitDict[c]])
                         aCluster = ROOT.sndCluster(c,1,hitlist,scifiDet,False)
                         clusScifi[index]=aCluster
                         index+=1
               cprev = c
          #print('number of clusters', ncl)
     return clusScifi
