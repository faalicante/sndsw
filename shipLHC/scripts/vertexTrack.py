import ROOT as r
#import rootUtils as ut
#import shipunit as u
from array import array
import ctypes
import time
import fedrarootlogon

r.gEDBDEBUGLEVEL = 0
r.gROOT.SetBatch(True)
h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=False)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=False)
parser.add_argument("--pdg", dest="nuPdg", help="nuPdg", required=True, type=int)
parser.add_argument("--offline", dest="OffMode", default=False, action="store_true")
parser.add_argument("--clusID", dest="ClusterID", required=False, default=0)
parser.add_argument("--procID", dest="ProcID", required=True, default=None, type=int)
options = parser.parse_args()

if options.OffMode:
     xroot_prefix = 'root:://eosuser.cern.ch/'
else:
     xroot_prefix = ''
pathGeofile = xroot_prefix+'/afs/cern.ch/work/f/falicant/public/matching'
#pathGeofile = '/home/fabio/Simulations/fedra'
pathSim = xroot_prefix+'/eos/user/a/aiulian/sim_fedra/numu_sim_activeemu_withcrisfiles_25_July_2022'
#pathSim = '/home/fabio/Simulations/fedra'
pathPlots = '/afs/cern.ch/work/f/falicant/public/vertex/histo_vtx/'
#pathPlots = '/home/fabio/Simulations/fedra'
pathText = '/afs/cern.ch/work/f/falicant/public/vertex/text_vtx/'
#pathText = '/home/fabio/Simulations/fedra'
nameGeofile = '/geofile_full.Genie-TGeant4.root'
nameSim = '/sndLHC.Genie-TGeant4.root'

brickID = ((options.ProcID//4)+1)*10+(options.ProcID%4)+1
nameVertex = '/b0000'+str(brickID)+'/vertextree.root'
#nameVertex = '/vertextree_31.root'


geoFile = pathGeofile+nameGeofile
simFile = pathSim+nameSim
vertexFile =pathSim+nameVertex

start_time = time.time()

fgeo = r.TFile.Open(geoFile)
#from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
#load geo dictionary
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
# -----Create geometry----------------------------------------------
import shipLHC_conf as sndDet_conf
import SndlhcGeo
geo = SndlhcGeo.GeoInterface(geoFile)
lsOfGlobals = r.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])

run = r.FairRunSim()
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()
nav = r.gGeoManager.GetCurrentNavigator()

scifiDet = r.gROOT.GetListOfGlobals().FindObject('Scifi')


'''
Wall_1               : z=  294.8486cm  dZ=    4.0747cm  [  290.7739  ,   298.9233] dx=   19.2996cm [  -45.9778   ,   -7.3787] dy=   19.2981cm [   15.1627   ,   53.7588]
     Brick_1         : z=  294.7286cm  dZ=    3.8947cm  [  290.8339   ,  298.6233] dx=    9.6498cm [  -26.6145    ,  -7.3149] dy=    9.6490cm [   15.2264    ,  34.5245] 
     Brick_2         : z=  294.7286cm  dZ=    3.8947cm  [  290.8339   ,  298.6233] dx=    9.6498cm [  -45.9141    , -26.6145] dy=    9.6490cm [   15.0989    ,  34.3970] 
     Brick_3         : z=  294.9686cm  dZ=    3.8947cm  [  291.0739   ,  298.8633] dx=    9.6498cm [  -26.7420    ,  -7.4424] dy=    9.6490cm [   34.5245    ,  53.8226]
     Brick_4         : z=  294.9686cm  dZ=    3.8947cm  [  291.0739   ,  298.8633] dx=    9.6498cm [  -46.0416    , -26.7420] dy=    9.6490cm [   34.3970    ,  53.6951]
Wall_2               : z=  307.8577cm  dZ=    4.0747cm  [  303.7830  ,   311.9324] dx=   19.2996cm [  -45.9668   ,   -7.3677] dy=   19.2981cm [   15.0110   ,   53.6071]
     Brick_1         : z=  307.7377cm  dZ=    3.8947cm  [  303.8430   ,  311.6324] dx=    9.6498cm [  -26.6035    ,  -7.3039] dy=    9.6490cm [   15.0747    ,  34.3728]
     Brick_2         : z=  307.7377cm  dZ=    3.8947cm  [  303.8430   ,  311.6324] dx=    9.6498cm [  -45.9031    , -26.6035] dy=    9.6490cm [   14.9472    ,  34.2453]
     Brick_3         : z=  307.9777cm  dZ=    3.8947cm  [  304.0830   ,  311.8724] dx=    9.6498cm [  -26.7310    ,  -7.4314] dy=    9.6490cm [   34.3728    ,  53.6709]
     Brick_4         : z=  307.9777cm  dZ=    3.8947cm  [  304.0830   ,  311.8724] dx=    9.6498cm [  -46.0306    , -26.7310] dy=    9.6490cm [   34.2453    ,  53.5434]
Wall_3               : z=  320.8470cm  dZ=    4.0747cm  [  316.7723  ,   324.9217] dx=   19.2996cm [  -45.9559   ,   -7.3567] dy=   19.2981cm [   14.8695   ,   53.4657]
     Brick_1         : z=  320.7270cm  dZ=    3.8947cm  [  316.8323   ,  324.6217] dx=    9.6498cm [  -26.5926    ,  -7.2930] dy=    9.6490cm [   14.9332    ,  34.2313]
     Brick_2         : z=  320.7270cm  dZ=    3.8947cm  [  316.8323   ,  324.6217] dx=    9.6498cm [  -45.8921    , -26.5926] dy=    9.6490cm [   14.8057    ,  34.1038]
     Brick_3         : z=  320.9670cm  dZ=    3.8947cm  [  317.0723   ,  324.8617] dx=    9.6498cm [  -26.7201    ,  -7.4205] dy=    9.6490cm [   34.2313    ,  53.5294]
     Brick_4         : z=  320.9670cm  dZ=    3.8947cm  [  317.0723   ,  324.8617] dx=    9.6498cm [  -46.0196    , -26.7201] dy=    9.6490cm [   34.1038    ,  53.4019]
Wall_4               : z=  333.8461cm  dZ=    4.0747cm  [  329.7714  ,   337.9208] dx=   19.2996cm [  -45.9549   ,   -7.3557] dy=   19.2981cm [   14.7178   ,   53.3140]
     Brick_1         : z=  333.7261cm  dZ=    3.8947cm  [  329.8314   ,  337.6208] dx=    9.6498cm [  -26.5916    ,  -7.2920] dy=    9.6490cm [   14.7816    ,  34.0797]
     Brick_2         : z=  333.7261cm  dZ=    3.8947cm  [  329.8314   ,  337.6208] dx=    9.6498cm [  -45.8911    , -26.5916] dy=    9.6490cm [   14.6541    ,  33.9522]
     Brick_3         : z=  333.9661cm  dZ=    3.8947cm  [  330.0714   ,  337.8608] dx=    9.6498cm [  -26.7191    ,  -7.4195] dy=    9.6490cm [   34.0797    ,  53.3778]
     Brick_4         : z=  333.9661cm  dZ=    3.8947cm  [  330.0714   ,  337.8608] dx=    9.6498cm [  -46.0186    , -26.7191] dy=    9.6490cm [   33.9522    ,  53.2503]
Wall_5               : z=  346.8452cm  dZ=    4.0747cm  [  342.7705  ,   350.9199] dx=   19.2996cm [  -45.9439   ,   -7.3447] dy=   19.2981cm [   14.5663   ,   53.1624]
     Brick_1         : z=  346.7252cm  dZ=    3.8947cm  [  342.8305   ,  350.6199] dx=    9.6498cm [  -26.5806    ,  -7.2810] dy=    9.6490cm [   14.6300    ,  33.9281]
     Brick_2         : z=  346.7252cm  dZ=    3.8947cm  [  342.8305   ,  350.6199] dx=    9.6498cm [  -45.8801    , -26.5806] dy=    9.6490cm [   14.5025    ,  33.8006]
     Brick_3         : z=  346.9652cm  dZ=    3.8947cm  [  343.0705   ,  350.8599] dx=    9.6498cm [  -26.7081    ,  -7.4085] dy=    9.6490cm [   33.9281    ,  53.2262]
     Brick_4         : z=  346.9652cm  dZ=    3.8947cm  [  343.0705   ,  350.8599] dx=    9.6498cm [  -46.0076    , -26.7081] dy=    9.6490cm [   33.8006    ,  53.0987]
'''
zConv = {'11':2986233, '12':2986233, '13':2988633, '14':2988633,
          '21':3116324, '22':3116324, '23':3118724, '24':3118724,
          '31':3246217, '32':3246217, '33':3248617, '34':3248617,
          '41':3376208, '42':3376208, '43':3378608, '44':3378608,
          '51':3506199, '52':3506199, '53':3508599, '54':3508599}

#brickID = ((ProcID//4)+1)*10+(ProcID%4)+1
#zConv[str(brickID)]


if options.nuPdg > 0:
	nu_pdg = options.nuPdg
	lep_pdg = options.nuPdg - 1
else:
	nu_pdg = options.nuPdg
	lep_pdg = options.nuPdg + 1

f1=r.TFile.Open(simFile)
eventTree = f1.cbmsim

dproc = r.EdbDataProc()
gAli = dproc.PVR()
scancond = r.EdbScanCond()
gAli.SetScanCond(scancond)

#setting parameters as in vertexing.C
vertexrec = r.EdbVertexRec()
vertexrec.SetPVRec(gAli)
vertexrec.eDZmax=3000.
vertexrec.eProbMin=0.0001
vertexrec.eImpMax=15.
vertexrec.eUseMom=False
vertexrec.eUseSegPar=True
vertexrec.eQualityMode=0

proc = r.EdbDataProc()

#reading vertices, building EdbVertex objects
print("start first loop to prepare vertices")

dproc.ReadVertexTree(vertexrec, vertexFile, "1")

vertices = gAli.eVTX
print("prepared vertices: ", len(vertices),"starting loop ")

#f2=r.TFile.Open(vertexFile)
#vertexTree = f2.vtx
#vertexTree.AddFriend(f2.brickinfo)

def GetCharge(pdgcode, pdgdb):
     charge = 0
     if pdgdb.GetParticle(pdgcode):
          charge = pdgdb.GetParticle(pdgcode).Charge()
     elif (pdgcode > 1e+8):
          charge = 1
     return charge

pdgdb = r.TDatabasePDG.Instance()

targetRangeX = [-50, 0]
targetRangeY = [10, 60]
targetRangeZ = [270, 370]
wallRangeX = [-46, -6]
wallRangeY = [14, 54]

single = False
histoFileName = pathPlots+str(options.ClusterID)+'_histoEml_'+str(brickID)+'.root'
#histoFileName = pathPlots+'/histoEml.root'

histoFile = r.TFile(histoFileName, 'recreate')
fitFileName = pathText+str(options.ClusterID)+'.'+str(nu_pdg)+'_fitEml_'+str(brickID)+'.txt'
#fitFileName = pathText+'/fitEml.txt'

statFileName = pathText+str(options.ClusterID)+'.'+str(nu_pdg)+'_statEml_'+str(brickID)+'.txt'
fitFile = open(fitFileName, "w+")
statFile = open(statFileName, "w+")
#fitFile.write('{0:<5}  {1:>10}  {2:>10}  {3:>10}  {4:>10}  {5:>10}  {6:>10}'.format('event','x_bar','x_rad','y_bar','y_rad','x_hits', 'y_hits')+'\n')
def addToFile(addPrint):
     fitFile.write(addPrint+'\n')

h_xy = r.TH2F('trk_xy', '2d emulsion map;x [cm]; y [cm]', 10000, targetRangeX[0], targetRangeX[1], 10000, targetRangeY[0], targetRangeY[1])

h_res_vtx_x = r.TH1F('res_vtx_x','x vtx residuals; x[um]',2000,-1000,1000)
h_res_vtx_y = r.TH1F('res_vtx_y','y vtx residuals; y[um]',2000,-1000,1000)
h_res_vtx_z = r.TH1F('res_vtx_z','z vtx residuals; y[um]',2000,-25000,50000)
h_res_vtx_d = r.TH1F('res_vtx_d','distance vtx residuals; d[um]',1000,0,1000)

h_res_x = r.TH1F('res_x','x residuals; x[cm]',200,-10,10)
h_res_y = r.TH1F('res_y','y residuals; y[cm]',200,-10,10)
h_res_d = r.TH1F('res_d','distance residuals; d[cm]',200,0,10)

A, B = r.TVector3(), r.TVector3()
ccCount=0
sameEvt = 0
sameWll = 0
count = 0
for vtx in vertices:
     if count%10==0: print(count/len(vertices))
     count+=1
     motherWall = -1
     moth = 0
     if vtx.Flag()==0 or vtx.Flag()==3:
          trackIDs = []
          eventIDs = []
          motherIDs = []
          eventDict = {}
          ntracks = vtx.N()
          for itrack in range(ntracks):
               track = vtx.GetTrack(itrack)
               trackIDs.append(track.MCTrack())
               eventIDs.append(track.MCEvt())
          #if eventIDs.count(eventIDs[0]) != len(eventIDs):
          #     sameEvt+=1
          #     continue
          for id in set(eventIDs):
               eventDict[id] = eventIDs.count(id)
          if max(eventDict.values()) > 0.5*ntracks:
               eventID = max(eventDict, key = eventDict.get)
          else:
               sameEvt+=1     
               continue
          vx = vtx.VX()
          vy = vtx.VY()
          vz = vtx.VZ()
          eventTree.GetEntry(eventID)
          incoming_nu = eventTree.MCTrack[0]
          outgoing_l = eventTree.MCTrack[1]
          for trackID in trackIDs:
               motherIDs.append(eventTree.MCTrack[trackID].GetMotherId())
          if motherIDs.count(0) < 0.5*len(motherIDs): continue
               #if eventTree.MCTrack[trackID].GetMotherId() == 0:
               #      moth = 1
          #if not moth: continue
          if abs(incoming_nu.GetPdgCode())!=nu_pdg or abs(outgoing_l.GetPdgCode())!=lep_pdg: continue
          mX = incoming_nu.GetStartX()
          mY = incoming_nu.GetStartY()
          mZ = incoming_nu.GetStartZ()
          r.gGeoManager.FindNode(mX, mY, mZ)
          node = r.gGeoManager.GetPath()
          for wall in range(5):
               if 'Wall_'+str(wall) in node:
                         motherWall = wall
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
          wallID = (options.ProcID//4)
          if motherWall != wallID:
               sameWll+=1
               continue
          ccCount+=1
          mX = mX * 1E+4 + 473000
          mY = mY * 1E+4 - 158000
          mZ = mZ * 1E+4 - zConv[str(brickID)]
          dist_vtx = r.TMath.Sqrt((vx-mX)**2+(vy-mY)**2)
          h_res_vtx_x.Fill(vx-mX)
          h_res_vtx_y.Fill(vy-mY)
          h_res_vtx_z.Fill(vz-mZ)
          h_res_vtx_d.Fill(dist_vtx)
          
          h_xy.Reset('ICES')
          for itrack in range(ntracks):
               track = vtx.GetTrack(itrack)
               nseg = track.N()
               seg = track.GetSegmentF(nseg-1)
               segPlate = seg.Plate()
               sx = seg.X()
               sy = seg.Y()
               if segPlate!=60:
                    sz = seg.Z()
                    tx = seg.TX()
                    ty = seg.TY()
                    sx += tx * (-sz)
                    sy += ty * (-sz)
               sx = (sx - 473000)/1E+4 
               sy = (sy + 158000)/1E+4 
               h_xy.Fill(sx, sy)
          barX = h_xy.GetMean(1)
          barY = h_xy.GetMean(2)
          radX = h_xy.GetStdDev(1)
          radY = h_xy.GetStdDev(2)
          dist = r.TMath.Sqrt((barX-nuX)**2+(barY-nuY)**2)
          h_res_x.Fill(barX-nuX)
          h_res_y.Fill(barY-nuY)
          h_res_d.Fill(dist)
          addToFile('{0:<5}  {1:>15}  {2:>15}  {3:>15}  {4:>15}  {5:>15}'.format(eventID, round(barX,2), round(radX,2), round(barY,2), round(radY,2), h_xy.GetEntries()))
print("Total number of event: {}".format(len(vertices)))
print("Selected number of event: {}".format(ccCount))
print("Different event: {}".format(sameEvt))
print("Different wall: {}".format(sameWll))
statFile.write("Total number of event: {}\n".format(len(vertices)))
statFile.write("Selected number of event: {}\n".format(ccCount))
statFile.write("Different event: {}\n".format(sameEvt))
statFile.write("Different wall: {}\n".format(sameWll))

histoFile.cd()
h_res_vtx_x.Write()
h_res_vtx_y.Write()
h_res_vtx_z.Write()
h_res_vtx_d.Write()
h_res_x.Write()
h_res_y.Write()
h_res_d.Write()

histoFile.Write()
histoFile.Close()
print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')
statFile.write('Done\n')
statFile.write('Elapsed time: '+str((time.time()-start_time)/60.)+' mins\n')
fitFile.close()
statFile.close()
