# read skimmed files, split them into different pt bins
import ROOT as r
import re
# if 'can' not in globals():
#     can = r.TCanvas('can', 'can', 600, 600)

# from subprocess import run, PIPE

# finlist = ["/home/tasheng/braa/CutSkim/BsMC.root", "/home/tasheng/braa/CutSkim/BsData.root"]
# finlist = ["BsMC.root", "BsData.root"]
pre = "/home/tasheng/braa/files/"

cat = "_newsample_presel_chi2all"
bpfinlist = [pre + post + cat + ".root" for post in ["BPData", "BPMC", "jpsinp"] ]
# bpptlist = [5, 7, 10, 15, 20, 50, 60]
# bpptlist = [5, 7, 10, 15, 20, 60]
# bpptlist = [5, 15, 60]
bpptlist = [5, 60]

bsfinlist = [pre + post + cat + ".root" for post in ["BsData", "BsMC"] ]
# bsptlist = [7, 10, 15, 20, 50]
# bsptlist = [7, 15, 50]
bsptlist = [7, 50]
# bsptlist = [10, 50]

chi2cut = 'Bchi2cl > 0.05'
dlscut = 'BsvpvDistRatio > 2'
cutlist = [chi2cut, dlscut]

for filename, treename, ptlist in zip([bpfinlist, bsfinlist],
                                      ['ntKp', 'ntphi'],
                                      [bpptlist, bsptlist]):
    for fstr in filename:
        for cut in cutlist:
            cutname = re.findall(r'\w+', cut)[0]
            dfnt = r.RDataFrame(treename, fstr)
            for pti, ptf in zip(ptlist[:-1], ptlist[1:]):
                dfnt.Filter(f'Bpt > {pti} && Bpt < {ptf}')\
                    .Filter(cut)\
                    .Snapshot(treename, fstr[:-5] + f'_{pti}_{ptf}_{cutname}.root')

# for fstr in bpfinlist:
#     fin = r.TFile(fstr)
#     dfnt = r.RDataFrame('ntKp', fin)
#     for pti, ptf in zip(bpptlist[:-1], bpptlist[1:]):
#         dfnt.Filter(f'Bpt > {pti} && Bpt < {ptf}').Snapshot(
#             'ntKp', fstr[:-5] + f'_{pti}_{ptf}.root')

# r.gStyle.SetOptStat(0)
# th2 = dfnt.Histo2D(( 'mod', ';Bchi2cl;DLS', 100, 0, 1, 100, 0, 80 ), 'Bchi2cl', 'BsvpvDistRatio')
# th2.SetTitle(f'Cov coef:{th2.GetCovariance():.2f}')
# th2.Draw('colorz')
# can.SaveAs('test.png')
