import ROOT


c = ROOT.TCanvas('residuals', 'residuals', 1200, 1200)
c.Divide(1, 2)
resX = ROOT.TH1D('resX', 'residuals x', 100, -20, 20)
resY = ROOT.TH1D('resY', 'residuals y', 100, -20, 20)

path = '/home/fabio/Simulations/numu_sim_25_July_2022/'
scifiFile = path+'scifi_fits.txt'
emulsionFile = path+'eml_fits.txt'
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
c.cd(1)
fitResX = resX.Fit('gaus','SQ','',-10,10)
fitResY = resY.Fit('gaus','SQ','',-10,10)
resMeanX = fitResX.Parameter(1)
resVarX = fitResX.Parameter(2)
resMeanY = fitResY.Parameter(1)
resVarY = fitResY.Parameter(2)
resX.Draw()
c.cd(2)
resY.Draw()
print(resMeanX, resVarX, resMeanY, resVarY)

c.Print('/home/fabio/Simulations/numu_sim_25_July_2022/eml_scifi_2.png')



