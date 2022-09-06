import ROOT
import rootUtils as ut
h = {}


path = '/home/fabio/Simulations/numu_sim_25_July_2022/'
scifiFile = path+'scifi_fits.txt'
emulsionFile = path+'eml_fits.txt'

def allEvents():
    c = ROOT.TCanvas('residuals', 'residuals')
    c.Divide(1, 3)
    c.SetCanvasSize(1000, 1500)
    c.SetWindowSize(1050, 1200)
    resX = ROOT.TH1D('resX', 'residuals x', 100, -20, 20)
    resY = ROOT.TH1D('resY', 'residuals y', 100, -20, 20)
    resD = ROOT.TH1D('resD', 'residuals distance', 100, 0, 20)

    with open(emulsionFile, 'r') as fe:
        linesEml = fe.readlines()
        for lineEml in linesEml:
            lineEml = lineEml.replace('\n', '')
            lineEml = lineEml.split()
            if lineEml[0] == 'event': continue
            emulsionEvent = int(lineEml[0])
            emulsionBarX = float(lineEml[1])
            emulsionBarY = float(lineEml[3])
            with open (scifiFile, 'r') as fs:
                linesScifi = fs.readlines()
                for lineScifi in linesScifi:
                    lineScifi = lineScifi.replace('\n', '')
                    lineScifi = lineScifi.split()
                    if lineScifi[0] == 'event': continue
                    scifiEvent = int(lineScifi[0])
                    if scifiEvent == emulsionEvent:
                        scifiBarX = float(lineScifi[1])
                        scifiBarY = float(lineScifi[3])
                        resX.Fill(emulsionBarX-scifiBarX)
                        resY.Fill(emulsionBarY-scifiBarY)
                        dist = ROOT.TMath.Sqrt((emulsionBarX-scifiBarX)**2+(emulsionBarY-scifiBarY)**2)
                        resD.Fill(dist)
    c.cd(1)
    fitResX = resX.Fit('gaus','SQ','',-10,10)
    fitResY = resY.Fit('gaus','SQ','',-10,10)
    resMeanX = fitResX.Parameter(1)
    resVarX = fitResX.Parameter(2)
    resMeanY = fitResY.Parameter(1)
    resVarY = fitResY.Parameter(2)
    eff0x = 0
    eff0y = 0
    eff1x = resX.GetEntries()
    eff1y = resY.GetEntries()
    for binX in range(resX.GetXaxis().GetNbins()):
        if resX.GetBinCenter(binX) < resMeanX - 3 * resVarX:
            eff0x += resX.GetBinContent(binX)
        elif resX.GetBinCenter(binX) > resMeanX + 3 * resVarX:
            eff0x += resX.GetBinContent(binX)
    for binY in range(resY.GetXaxis().GetNbins()):
        if resY.GetBinCenter(binY) < resMeanY - 3 * resVarY:
            eff0y += resY.GetBinContent(binY)
        elif resY.GetBinCenter(binY) > resMeanY + 3 * resVarY:
            eff0y += resY.GetBinContent(binY)
    effx = 1-eff0x/eff1x
    effy = 1-eff0y/eff1y
    print('efficiency x:', effx, eff0x, eff1x)
    print('efficiency y:', effy, eff0y, eff1y)
    resX.Draw()
    c.cd(2)
    resY.Draw()
    c.cd(3)
    resD.Draw()
    print(resMeanX, resVarX, resMeanY, resVarY)

    c.Print('/home/fabio/Simulations/numu_sim_25_July_2022/eml_scifi_2.png')

def blockEvents(Nint = 200):
    #c = ROOT.TCanvas('resid', 'residuals', 1280, 720)
    #distT = ROOT.TH1D('distT', 'true event distance', 10, 0, 1)
    #distF = ROOT.TH1D('distF', 'false event distance', 10, 0, 1)
    tCount = 0
    fCount = 0
    block = 0
    average = 0
    with open(emulsionFile, 'r') as fe:
        N1 = 0
        linesEml = fe.readlines()
        for lineEml in linesEml:
            lineEml = lineEml.replace('\n', '')
            lineEml = lineEml.split()
            if lineEml[0] == 'event': continue
            N1+=1
            emulsionEvent = int(lineEml[0])
            emulsionBarX = float(lineEml[1])
            emulsionBarY = float(lineEml[3])
            #print('eml', emulsionEvent)
            minDist = 100
            with open (scifiFile, 'r') as fs:
                N2 = 0
                linesScifi = fs.readlines()
                if Nint*(block+1) > len(linesScifi): break
                for x in range(Nint):
                    y = x + Nint * block
                    linesScifi[y] = linesScifi[y].replace('\n', '')
                    linesScifi[y] = linesScifi[y].split()
                    N2+=1
                    scifiBarX = float(linesScifi[y][1])
                    scifiBarY = float(linesScifi[y][3])
                    dist = ROOT.TMath.Sqrt((emulsionBarX-scifiBarX)**2+(emulsionBarY-scifiBarY)**2)
                    #print('scifi', int(linesScifi[y][0]))

                    #print(dist)
                    if dist < minDist:
                        minDist = dist
                        scifiEvent = int(linesScifi[y][0])
                if scifiEvent == emulsionEvent:
                    #distT.Fill(minDist)
                    tCount+=1
                else:
                    #distF.Fill(minDist)
                    fCount+=1
            if N1 == Nint: 
                N1 = 0
                block+=1
                #tCount = distT.GetEntries()
                #fCount = distF.GetEntries()
                print('block {} true {} false {}'.format(block, tCount, fCount))
                print('block {} matching efficiency {}'.format(block, tCount/(tCount+fCount)))
                average += tCount/(tCount+fCount)
                tCount = 0
                fCount = 0
    average = average/block
    print('average efficiency', average)      


    
    #distT.SetMaximum(100)
    #distT.SetLineColor(4)
    #distF.SetLineColor(2)
    #distT.Draw()
    #distF.Draw('same')
    #c.Print('/home/fabio/Scrivania/match.png')

def hitCorr():
    hits = []
    tracks = []
    with open(emulsionFile, 'r') as fe:
        linesEml = fe.readlines()
        for lineEml in linesEml:
            lineEml = lineEml.replace('\n', '')
            lineEml = lineEml.split()
            if lineEml[0] == 'event': continue
            tracks.append(float(lineEml[5]))
    with open (scifiFile, 'r') as fs:
        linesScifi = fs.readlines()
        for lineScifi in linesScifi:
            lineScifi = lineScifi.replace('\n', '')
            lineScifi = lineScifi.split()
            if lineScifi[0] == 'event': continue
            hits.append(float(lineScifi[5])+float(lineScifi[6]))
    minTrack = min(tracks)-1
    maxTrack = max(tracks)+1
    minHit = min(hits)-1
    maxHit = max(hits)+1
    print(minTrack, maxTrack, minHit, maxHit)
    entries = len(tracks)
    if len(tracks) != len(hits):
        print('Error: different size')
    ut.bookCanvas(h, 'corr', 'correlation', 1280, 720)
    ut.bookHist(h, 'hits', 'hits vs tracks;#tracks;#hits', 100, minTrack, maxTrack, 100, minHit, maxHit)
    for i in range(entries):
        h['hits'].Fill(tracks[i], hits[i])
    h['hits'].SetMarkerStyle(20)
    h['hits'].SetMarkerColor(1)
    h['hits'].SetMarkerSize(0.4)
    h['hits'].Draw()
    h['corr'].Print('/home/fabio/cernbox/softphys/corr_all.png')