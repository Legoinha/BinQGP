# read skimmed files, split them into different pt bins
import ROOT as r

# from subprocess import run, PIPE

# finlist = ["/home/tasheng/braa/CutSkim/BsMC.root", "/home/tasheng/braa/CutSkim/BsData.root"]
# finlist = ["BsMC.root", "BsData.root"]
pre = "/eos/user/h/hmarques/work/data/"

bpfinlist = [pre + post + ".root" for post in ["BPData_nom", "BPMC_nom"] ]
bpptlist = [5, 7, 10, 15, 20, 30, 50, 60]

bsfinlist = [pre + post + ".root" for post in ["BsData_nom", "BsMC_nom"] ]
bsptlist = [7, 10, 15, 20, 50]

for fstr in bsfinlist:
    fin = r.TFile(fstr)
    dfnt = r.RDataFrame('ntphi', fin)
    for pti, ptf in zip(bsptlist[:-1], bsptlist[1:]):
        dfnt.Filter(f'Bpt > {pti} && Bpt < {ptf}').Snapshot('ntphi', fstr[:-5] + f'_{pti}_{ptf}.root')

for fstr in bpfinlist:
    fin = r.TFile(fstr)
    dfnt = r.RDataFrame('ntKp', fin)
    for pti, ptf in zip(bpptlist[:-1], bpptlist[1:]):
        dfnt.Filter(f'Bpt > {pti} && Bpt < {ptf}').Snapshot('ntKp', fstr[:-5] + f'_{pti}_{ptf}.root')
