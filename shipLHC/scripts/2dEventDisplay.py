import ROOT,os
import rootUtils as ut
from array import array
import shipunit as u
import SndlhcMuonReco
A,B = ROOT.TVector3(),ROOT.TVector3()
freq      =  160.316E6

h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=False)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input file MC",default="",required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)

parser.add_argument("-H", "--houghTransform", dest="houghTransform", help="do not use hough transform for track reco", action='store_false',default=True)
parser.add_argument("-t", "--tolerance", dest="tolerance",  type=float, help="How far away from Hough line hits assigned to the muon can be. In cm.", default=0.)
parser.add_argument("--hits_to_fit", dest = "hits_to_fit", type=str, help="Which detectors to use in the fit, in the format: vesfusds, where [ve] is veto, [sf] is Scifi, [us] is Upstream muon filter, and [ds] is downstream muon filter", default = "sfusds")
parser.add_argument("--hits_for_triplet", dest = "hits_for_triplet", type=str, help="Which detectors to use for the triplet condition. In the same format as --hits_to_fit", default = "ds")

options = parser.parse_args()

trans2local = False

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geoFile)

lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])
lsOfGlobals.Add(geo.modules['MuFilter'])

detSize = {}
si = geo.snd_geo.Scifi
detSize[0] =[si.channel_width, si.channel_width, si.scifimat_z ]
mi = geo.snd_geo.MuFilter
detSize[1] =[mi.VetoBarX/2,                   mi.VetoBarY/2,            mi.VetoBarZ/2]
detSize[2] =[mi.UpstreamBarX/2,           mi.UpstreamBarY/2,    mi.UpstreamBarZ/2]
detSize[3] =[mi.DownstreamBarX_ver/2,mi.DownstreamBarY/2,mi.DownstreamBarZ/2]

mc = False

run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()

if options.inputFile=="":
  f=ROOT.TFile.Open(options.path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
else:
  f=ROOT.TFile.Open(options.path+options.inputFile)

if f.FindKey('cbmsim'):
        eventTree = f.cbmsim
        runId = 'sim'
        if eventTree.GetBranch('ScifiPoint'): mc = True
else:   
        eventTree = f.rawConv
        ioman.SetTreeName('rawConv')
        #runId = eventTree.EventHeader.GetRunId()

outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(f)
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

if options.houghTransform:
  muon_reco_task = SndlhcMuonReco.MuonReco()
  run.AddTask(muon_reco_task)
else:
  import SndlhcTracking
  trackTask = SndlhcTracking.Tracking() 
  trackTask.SetName('simpleTracking')
  run.AddTask(trackTask)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

run.Init()
eventTree = ioman.GetInTree()
# backward compatbility for early converted events
eventTree.GetEvent(0)
if eventTree.GetBranch('Digi_MuFilterHit'): eventTree.Digi_MuFilterHits = eventTree.Digi_MuFilterHit

if options.houghTransform:
# prepare track reco with hough transform
  muon_reco_task.SetTolerance(options.tolerance)
  muon_reco_task.SetHitsToFit(options.hits_to_fit)
  muon_reco_task.SetHitsForTriplet(options.hits_for_triplet)

nav = ROOT.gGeoManager.GetCurrentNavigator()

Nlimit = 4
onlyScifi = False
def goodEvent(event):
# can be replaced by any user selection
           stations = {'Scifi':{},'Mufi':{}}
           if event.Digi_ScifiHits.GetEntries()>25: return False
           for d in event.Digi_ScifiHits:
               stations['Scifi'][d.GetDetectorID()//1000000] = 1
           for d in event.Digi_MuFilterHits:
               plane = d.GetDetectorID()//1000
               stations['Mufi'][plane] = 1
           totalN = len(stations['Mufi'])+len(stations['Scifi'])
           if len(stations['Scifi'])>4 and len(stations['Mufi'])>6: return True
           else: False
           if onlyScifi and len(stations['Scifi'])>Nlimit: return True
           elif not onlyScifi  and totalN >  Nlimit: return True
           else: return False

def loopEvents(start=0,save=False,goodEvents=False,withTrack=-1,nTracks=0,minSipmMult=1, Setup='',verbose=0):
 if 'simpleDisplay' not in h: ut.bookCanvas(h,key='simpleDisplay',title='simple event display',nx=1200,ny=1600,cx=1,cy=2)
 h['simpleDisplay'].cd(1)
 zStart = 250. # TI18 coordinate system
 if Setup == 'H6': zStart = 60.
 if Setup == 'TP': zStart = -50. # old coordinate system with origin in middle of target
 if 'xz' in h: 
        h.pop('xz').Delete()
        h.pop('yz').Delete()
 ut.bookHist(h,'xz','; z [cm]; x [cm]',500,zStart,zStart+350.,100,-100.,10.)
 ut.bookHist(h,'yz','; z [cm]; y [cm]',500,zStart,zStart+350.,100,-30.,80.)

 proj = {1:'xz',2:'yz'}
 h['xz'].SetStats(0)
 h['yz'].SetStats(0)

 N = -1
 Tprev = -1
 A,B = ROOT.TVector3(),ROOT.TVector3()
 ptext={0:'   Y projection',1:'   X projection'}
 text = ROOT.TLatex()
 event = eventTree
 OT = sink.GetOutTree()
 if withTrack==0: OT = eventTree
 for N in range(start, event.GetEntries()):
    rc = event.GetEvent(N)
    if goodEvents and not goodEvent(event): continue
    if not withTrack<-1:
       if options.houghTransform:
          OT.Reco_MuonTracks.Delete()
          rc = source.GetInTree().GetEvent(N)
          muon_reco_task.Exec(0)
          ntracks = OT.Reco_MuonTracks.GetEntries()
          uniqueTracks = cleanTracks()
          if len(uniqueTracks)<nTracks: continue
       else:
          if withTrack==2:  trackTask.ExecuteTask("Scifi")
          elif withTrack==3:  trackTask.ExecuteTask("DS")
       ntracks = len(OT.Reco_MuonTracks)
       if ntracks<nTracks: continue

       if verbose>0:
          for aTrack in OT.Reco_MuonTracks:
             mom    = aTrack.getFittedState().getMom()
             pos      = aTrack.getFittedState().getPos()
             mom.Print()
             pos.Print()
    T,dT = 0,0
    if event.FindBranch("EventHeader"):
       T = event.EventHeader.GetEventTime()
       runId = eventTree.EventHeader.GetRunId()
       if Tprev >0: dT = T-Tprev
       Tprev = T
    if ntracks > 0: print('number of tracks: ', ntracks)

    digis = []
    if event.FindBranch("Digi_ScifiHits"): digis.append(event.Digi_ScifiHits)
    if event.FindBranch("Digi_MuFilterHits"): digis.append(event.Digi_MuFilterHits)
    if event.FindBranch("Digi_MuFilterHit"): digis.append(event.Digi_MuFilterHit)
    empty = True
    for x in digis:
       if x.GetEntries()>0:
         empty = False
         print( "event -> %i"%N)
    if empty: continue
    h['hitCollectionX']= {'Scifi':[0,ROOT.TGraphErrors()],'DS':[0,ROOT.TGraphErrors()]}
    h['hitCollectionY']= {'Veto':[0,ROOT.TGraphErrors()],'Scifi':[0,ROOT.TGraphErrors()],'US':[0,ROOT.TGraphErrors()],'DS':[0,ROOT.TGraphErrors()]}
    h['firedChannelsX']= {'Scifi':[0,0,0],'DS':[0,0,0]}
    h['firedChannelsY']= {'Veto':[0,0,0,0],'Scifi':[0,0,0],'US':[0,0,0,0],'DS':[0,0,0,0]}
    systems = {1:'Veto',2:'US',3:'DS',0:'Scifi'}
    for collection in ['hitCollectionX','hitCollectionY']:
       for c in h[collection]:
          rc=h[collection][c][1].SetName(c)
          rc=h[collection][c][1].Set(0)

    dTs = "%5.2Fns"%(dT/freq*1E9)
    # find detector which triggered
    minT = firstTimeStamp(event)
    dTs+= "    " + str(minT[1].GetDetectorID())
    for p in proj:
       rc = h[ 'simpleDisplay'].cd(p)
       h[proj[p]].Draw('b')
    emptyNodes()
    drawDetectors()
    for D in digis:
      for digi in D:
         detID = digi.GetDetectorID()
         sipmMult = 1
         if digi.GetName()  == 'MuFilterHit':
            system = digi.GetSystem()
            geo.modules['MuFilter'].GetPosition(detID,A,B)
            sipmMult = len(digi.GetAllSignals())
            if sipmMult<minSipmMult and (system==1 or system==2): continue
            curPath = nav.GetPath()
            tmp = curPath.rfind('/')
            nav.cd(curPath[:tmp])
         else:
            geo.modules['Scifi'].GetSiPMPosition(detID,A,B)
            curPath = nav.GetPath()
            tmp = curPath.rfind('/')
            nav.cd(curPath[:tmp])
            system = 0
         globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
         if trans2local:   nav.MasterToLocal(globA,locA)
         Z = A[2]
         if digi.isVertical():
                   collection = 'hitCollectionX'
                   Y = locA[0]
                   sY = detSize[system][0]
         else:                         
                   collection = 'hitCollectionY'
                   Y = locA[1]
                   sY = detSize[system][1]
         c = h[collection][systems[system]]
         rc = c[1].SetPoint(c[0],Z, Y)
         rc = c[1].SetPointError(c[0],detSize[system][2],sY)
         c[0]+=1 
         fillNode(curPath)
         if digi.isVertical():  F = 'firedChannelsX'
         else:                     F = 'firedChannelsY'
         ns = max(1,digi.GetnSides())
         for side in range(ns):
             for m in  range(digi.GetnSiPMs()):
                   qdc = digi.GetSignal(m+side*digi.GetnSiPMs())
                   if qdc < 0 and qdc > -900:  h[F][systems[system]][1]+=1
                   elif not qdc<0:   
                       h[F][systems[system]][0]+=1
                       #h[F][systems[system]][2+side]+=qdc
    h['hitCollectionY']['Scifi'][1].SetMarkerColor(ROOT.kBlue+2)
    h['hitCollectionX']['Scifi'][1].SetMarkerColor(ROOT.kBlue+2)
    k = 1
    for collection in ['hitCollectionX','hitCollectionY']:
       h[ 'simpleDisplay'].cd(k)
       drawInfo(h[ 'simpleDisplay'], k, runId, N, T)
       k+=1
       for c in h[collection]:
          F = collection.replace('hitCollection','firedChannels')
          pj = collection.split('ion')[1]
          if pj =="X" or c=="Scifi":
              print( "%1s %5s %3i  +:%3i -:%3i qdc :%5.1F"%(pj,c,h[collection][c][1].GetN(),h[F][c][0],h[F][c][1],h[F][c][2]))
          else:
              print( "%1s %5s %3i  +:%3i -:%3i qdcL:%5.1F qdcR:%5.1F"%(pj,c,h[collection][c][1].GetN(),h[F][c][0],h[F][c][1],h[F][c][2],h[F][c][3]))
          if h[collection][c][1].GetN()<1: continue
          if c=='Scifi':
            h[collection][c][1].SetMarkerStyle(20)
            h[collection][c][1].SetMarkerSize(1.5)
            rc=h[collection][c][1].Draw('sameP')
            h['display:'+c]=h[collection][c][1]
    h[ 'simpleDisplay'].Update()
    if withTrack == 2: addTrack(OT,True) #withTrack=2 scifi, =3 DS
    elif not withTrack<0:  addTrack(OT)
    if verbose>0: dumpChannels()
    if save: h['simpleDisplay'].Print('{:0>2d}-event_{:04d}'.format(runId,N)+'.png')
    rc = input("hit return for next event or q for quit: ")
    if rc=='q': break
 if save: os.system("convert -delay 60 -loop 0 event*.png animated.gif")

def addTrack(OT,scifi=False):
   xax = h['xz'].GetXaxis()
   nTrack = 0
   for   aTrack in OT.Reco_MuonTracks:
      trackColor = ROOT.kRed
      if aTrack.GetUniqueID()==1: trackColor = ROOT.kBlue+2
      if aTrack.GetUniqueID()==3: trackColor = ROOT.kBlack
      S = aTrack.getFitStatus()
      if not S.isFitConverged() and scifi:
         print('not converge')
         continue
      for p in [0,1]:
          h['aLine'+str(nTrack*10+p)] = ROOT.TGraph()

      zEx = xax.GetBinCenter(1)
      mom    = aTrack.getFittedState().getMom()
      pos      = aTrack.getFittedState().getPos()
      lam      = (zEx-pos.z())/mom.z()
      Ex        = [pos.x()+lam*mom.x(),pos.y()+lam*mom.y()]
      for p in [0,1]:   h['aLine'+str(nTrack*10+p)].SetPoint(0,zEx,Ex[p])

      for i in range(aTrack.getNumPointsWithMeasurement()):
         state = aTrack.getFittedState(i)
         pos    = state.getPos()
         for p in [0,1]:
             h['aLine'+str(nTrack*10+p)].SetPoint(i+1,pos[2],pos[p])

      zEx = xax.GetBinCenter(xax.GetLast())
      mom    = aTrack.getFittedState().getMom()
      pos      = aTrack.getFittedState().getPos()
      lam      = (zEx-pos.z())/mom.z()
      Ex        = [pos.x()+lam*mom.x(),pos.y()+lam*mom.y()]
      for p in [0,1]:   h['aLine'+str(nTrack*10+p)].SetPoint(i+2,zEx,Ex[p])

      for p in [0,1]:
             tc = h[ 'simpleDisplay'].cd(p+1)
             h['aLine'+str(nTrack*10+p)].SetLineColor(trackColor)
             h['aLine'+str(nTrack*10+p)].SetLineWidth(2)
             h['aLine'+str(nTrack*10+p)].Draw('same')
             tc.Update()
             h[ 'simpleDisplay'].Update()
      nTrack+=1

def drawDetectors():
   nodes = {'volMuFilter_1/volFeBlockEnd_1':ROOT.kGreen-3}
   for i in range(2):
      nodes['volVeto_1/volVetoPlane_{}_{}'.format(i, i)]=ROOT.kRed
      for j in range(7):
         nodes['volVeto_1/volVetoPlane_{}_{}/volVetoBar_1{}{:0>3d}'.format(i, i, i, j)]=ROOT.kRed
      nodes['volVeto_1/subVetoBox_{}'.format(i)]=ROOT.kGray+1
   for i in range(4):
      nodes['volMuFilter_1/volMuDownstreamDet_{}_{}'.format(i, i+7)]=ROOT.kBlue+2
      for j in range(60):
         nodes['volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_ver_3{}{:0>3d}'.format(i, i+7, i, j+60)]=ROOT.kBlue+2
         if i < 3:
            nodes['volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_hor_3{}{:0>3d}'.format(i, i+7, i, j)]=ROOT.kBlue+2
   for i in range(4):
      nodes['volMuFilter_1/subDSBox_{}'.format(i+7)]=ROOT.kGray
   for i in range(5):
      nodes['volTarget_1/ScifiVolume{}_{}000000'.format(i+1, i+1)]=ROOT.kBlue+2
      nodes['volTarget_1/volWallborder_{}'.format(i)]=ROOT.kGray
      nodes['volMuFilter_1/subUSBox_{}'.format(i+2)]=ROOT.kGray
      nodes['volMuFilter_1/volMuUpstreamDet_{}_{}'.format(i, i+2)]=ROOT.kBlue+2
      for j in range(10):
         nodes['volMuFilter_1/volMuUpstreamDet_{}_{}/volMuUpstreamBar_2{}00{}'.format(i, i+2, i, j)]=ROOT.kBlue+2
      nodes['volMuFilter_1/volFeBlock_{}'.format(i)]=ROOT.kGreen-3
   for i in range(7,10):
      nodes['volMuFilter_1/volFeBlock_{}'.format(i)]=ROOT.kGreen-3
   passNodes = {'Block', 'Wall'}
   xNodes = {'UpstreamBar', 'VetoBar', 'hor'}
   proj = {'X':0,'Y':1}
   for node_ in nodes:
      node = '/cave_1/Detector_0/'+node_
      for p in proj:
         if node+p not in h:
            nav.cd(node)
            N = nav.GetCurrentNode()
            S = N.GetVolume().GetShape()
            dx,dy,dz = S.GetDX(),S.GetDY(),S.GetDZ()
            ox,oy,oz = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
            P = {}
            M = {}
            if p=='X' and not any(xNode in node for xNode in xNodes):
               P['LeftBottom'] = array('d',[-dx+ox,oy,-dz+oz])
               P['LeftTop'] = array('d',[dx+ox,oy,-dz+oz])
               P['RightBottom'] = array('d',[-dx+ox,oy,dz+oz])
               P['RightTop'] = array('d',[dx+ox,oy,dz+oz])
            elif p=='Y' and 'ver' not in node:
               P['LeftBottom'] = array('d',[ox,-dy+oy,-dz+oz])
               P['LeftTop'] = array('d',[ox,dy+oy,-dz+oz])
               P['RightBottom'] = array('d',[ox,-dy+oy,dz+oz])
               P['RightTop'] = array('d',[ox,dy+oy,dz+oz])
            else: continue
            for C in P:
               M[C] = array('d',[0,0,0])
               nav.LocalToMaster(P[C],M[C])
            h[node+p] = ROOT.TPolyLine()
            X = h[node+p]
            c = proj[p]
            X.SetPoint(0,M['LeftBottom'][2],M['LeftBottom'][c])
            X.SetPoint(1,M['LeftTop'][2],M['LeftTop'][c])
            X.SetPoint(2,M['RightTop'][2],M['RightTop'][c])
            X.SetPoint(3,M['RightBottom'][2],M['RightBottom'][c])
            X.SetPoint(4,M['LeftBottom'][2],M['LeftBottom'][c])
            X.SetLineColor(nodes[node_])
            X.SetLineWidth(0)
            h[ 'simpleDisplay'].cd(c+1)
            if any(passNode in node for passNode in passNodes):
               X.SetFillColorAlpha(nodes[node_], 0.5)
               X.Draw('f&&same')
            X.Draw('same')
         else:
            X = h[node+p]
            c = proj[p]
            h[ 'simpleDisplay'].cd(c+1)
            if any(passNode in node for passNode in passNodes):
               X.Draw('f&&same') 
            X.Draw('same')

def dumpVeto():
    muHits = {10:[],11:[]}
    for aHit in eventTree.Digi_MuFilterHits:
         if not aHit.isValid(): continue
         s = aHit.GetDetectorID()//10000
         if s>1: continue
         p = (aHit.GetDetectorID()//1000)%10
         bar = (aHit.GetDetectorID()%1000)%60
         plane = s*10+p
         muHits[plane].append(aHit)
    for plane in [10,11]:
        for aHit in muHits[plane]:
          S =aHit.GetAllSignals()
          txt = ""
          for x in S:
              if x[1]>0: txt+=str(x[1])+" "
          print(plane, (aHit.GetDetectorID()%1000)%60, txt)

def cleanTracks():
    OT = sink.GetOutTree()
    listOfDetIDs = {}
    n = 0
    for aTrack in OT.Reco_MuonTracks:
        listOfDetIDs[n] = []
        for i in range(aTrack.getNumPointsWithMeasurement()):
           M =  aTrack.getPointWithMeasurement(i)
           R =  M.getRawMeasurement()
           listOfDetIDs[n].append(R.getDetId())
           if R.getDetId()>0: listOfDetIDs[n].append(R.getDetId()-1)
           listOfDetIDs[n].append(R.getDetId()+1)
        n+=1
    uniqueTracks = []
    for n1 in range( len(listOfDetIDs) ):
       unique = True
       for n2 in range( len(listOfDetIDs) ):
             if n1==n2: continue
             I = set(listOfDetIDs[n1]).intersection(listOfDetIDs[n2])
             if len(I)>0:  unique = False
       if unique: uniqueTracks.append(n1)
    if len(uniqueTracks)>1: 
         for n1 in range( len(listOfDetIDs) ): print(listOfDetIDs[n1])
    return uniqueTracks

def mufiNoise():
  for s in range(1,4): 
    ut.bookHist(h,'mult'+str(s),'hit mult for system '+str(s),100,-0.5,99.5)
    ut.bookHist(h,'multb'+str(s),'hit mult per bar for system '+str(s),20,-0.5,19.5)
    ut.bookHist(h,'res'+str(s),'residual system '+str(s),20,-10.,10.)
  OT = sink.GetOutTree()
  N=0
  for event in eventTree:
       N+=1
       if N%1000==0: print(N)
       OT.Reco_MuonTracks.Delete()
       rc = trackTask.ExecuteTask("Scifi")
       for aTrack in OT.Reco_MuonTracks:
           mom    = aTrack.getFittedState().getMom()
           pos      = aTrack.getFittedState().getPos()
           if not aTrack.getFitStatus().isFitConverged(): continue
           mult = {1:0,2:0,3:0}
           for aHit in eventTree.Digi_MuFilterHits:
              if not aHit.isValid(): continue
              s = aHit.GetDetectorID()//10000
              S = aHit.GetAllSignals()
              rc = h['multb'+str(s)].Fill(len(S))
              mult[s]+=len(S)
              if s==2 or s==1:
                 geo.modules['MuFilter'].GetPosition(aHit.GetDetectorID(),A,B)
                 y = (A[1]+B[1])/2.
                 zEx = (A[2]+B[2])/2.
                 lam      = (zEx-pos.z())/mom.z()
                 Ey        = pos.y()+lam*mom.y()
                 rc = h['res'+str(s)].Fill(Ey-y)
           for s in mult: rc = h['mult'+str(s)].Fill(mult[s])
  ut.bookCanvas(h,'noise','',1200,1200,2,3)
  for s in range(1,4):
   tc = h['noise'].cd(s*2-1)
   tc.SetLogy(1)
   h['mult'+str(s)].Draw()
   h['noise'].cd(s*2)
   h['multb'+str(s)].Draw()
  ut.bookCanvas(h,'res','',600,1200,1,3)
  for s in range(1,4):
   tc = h['res'].cd(s)
   h['res'+str(s)].Draw()

def firstTimeStamp(event):
        tmin = [1E9,'']
        digis = [event.Digi_MuFilterHits,event.Digi_ScifiHits]
        for digi in event.Digi_ScifiHits:
               dt = digi.GetTime()
               if dt<tmin[0]:
                    tmin[0]=dt
                    tmin[1]=digi
        for digi in event.Digi_MuFilterHits:
           for t in digi.GetAllTimes():
               dt = t.second
               if dt<tmin[0]:
                    tmin[0]=dt
                    tmin[1]=digi
        return tmin

def dumpChannels(D='Digi_MuFilterHits'):
     X = eval("eventTree."+D)
     text = {}
     for aHit in X:
         side = 'L'
         txt = "%8i"%(aHit.GetDetectorID())
         for k in range(aHit.GetnSiPMs()*aHit.GetnSides()):
              qdc = aHit.GetSignal(k)
              if qdc < -900: continue
              i = k
              if not k<aHit.GetnSiPMs():
                   i = k-aHit.GetnSiPMs()
                   if side == 'L':
                        txt += " | "
                        side = 'R'
              txt+= "  %2i:%4.1F  "%(i,qdc)
         text[aHit.GetDetectorID()] = txt
     keys = list(text.keys())
     keys.sort()
     for k in keys: print(text[k])

def fillNode(node):
   xNodes = {'UpstreamBar', 'VetoBar', 'hor'}
   proj = {'X':0,'Y':1}
   color = ROOT.kBlack
   thick = 10
   for p in proj:
      if node+p in h:
         X = h[node+p]
         if 'Veto' in node:
            color = ROOT.kRed+1
         if 'Downstream' in node:
            thick = 5
         c = proj[p]
         h[ 'simpleDisplay'].cd(c+1)
         X.SetFillColor(color)
         X.SetLineColor(color)
         X.SetLineWidth(thick)
         X.Draw('f&&same')
         X.Draw('same')   

def emptyNodes():
   nodes = {}
   for i in range(2):
      for j in range(7):
         nodes['volVeto_1/volVetoPlane_{}_{}/volVetoBar_1{}{:0>3d}'.format(i, i, i, j)]=ROOT.kRed
   for i in range(3):
      for j in range(60):
         nodes['volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_hor_3{}{:0>3d}'.format(i, i+7, i, j)]=ROOT.kBlue+2
         nodes['volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_ver_3{}{:0>3d}'.format(i, i+7, i, j+60)]=ROOT.kBlue+2
   for j in range(60):
      nodes['volMuFilter_1/volMuDownstreamDet_3_10/volMuDownstreamBar_ver_33{:0>3d}'.format(j+60)]=ROOT.kBlue+2
   for i in range(5):
      for j in range(10):
         nodes['volMuFilter_1/volMuUpstreamDet_{}_{}/volMuUpstreamBar_2{}00{}'.format(i, i+2, i, j)]=ROOT.kBlue+2
   proj = {'X':0,'Y':1}
   for node_ in nodes:
      node = '/cave_1/Detector_0/'+node_
      for p in proj:
         try:
            X = h[node+p]
            if X.GetFillColor()!=19:
               c = proj[p]
               h[ 'simpleDisplay'].cd(c+1)
               X.SetFillColorAlpha(19, 0)
               X.SetLineColor(nodes[node_])
               X.SetLineWidth(1)
               X.Draw('f&&same')
               X.Draw('same')
         except:
            notFilled = 1

def drawInfo(pad, k, run, event, time):
   drawLogo = True
   drawText = True
   if drawLogo:
      padLogo = ROOT.TPad("logo","logo",0.1,0.1,0.2,0.3)
      padLogo.SetFillStyle(4000)
      padLogo.SetFillColorAlpha(0, 0)
      padLogo.Draw()
      logo = ROOT.TImage.Open('/home/fabio/Immagini/Large__SND_Logo_black_cut.png')
      logo.SetConstRatio(True)
      logo.DrawText(0, 0, 'SND', 98)
      padLogo.cd()
      logo.Draw()
      pad.cd(k)

   if drawText:
      padText = ROOT.TPad("info","info",0.19,0.1,0.36,0.32)
      padText.SetFillStyle(4000)
      padText.Draw()
      padText.cd()
      textInfo = ROOT.TLatex()
      textInfo.SetTextAlign(11)
      textInfo.SetTextFont(42)
      textInfo.SetTextSize(.15)
      textInfo.DrawLatex(0, 0.6, 'SND@LHC Experiment, CERN')
      textInfo.DrawLatex(0, 0.4, 'Run / Event: '+str(run)+' / '+str(event))
      textInfo.DrawLatex(0, 0.2, 'Time Stamp: {} a.u.'.format(time))
      pad.cd(k)
