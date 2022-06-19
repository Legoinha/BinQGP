#! /bin/python
import ROOT as r
import shutil
ptlist = [5, 7, 10, 15, 20, 50]

outFiles = ['BPw.root', 'Bsw.root']

fout = r.TFile(outFiles[0], 'recreate')
for pti, ptf in zip(ptlist[:-1], ptlist[1:]):
    fin = r.TFile(f'./results/Bu/{pti}_{ptf}/mc_validation_plots/weights/weights.root')
    h = fin.Get(f'weights_BDT_pt_{pti}_{ptf}')
    fout.cd()
    h.Write()
fout.Close()

ptlist = [7, 10, 15, 20]
fout = r.TFile(outFiles[1], 'recreate')
for pti, ptf in zip(ptlist[:-1], ptlist[1:]):
    fin = r.TFile(f'./results/Bs/{pti}_{ptf}/mc_validation_plots/weights/weights.root')
    h = fin.Get(f'weights_BDT_pt_{pti}_{ptf}')
    fout.cd()
    h.Write()
fout.Close()

shutil.copy2(outFiles[0], '../gd/braa_nohbhe_trk/BP/EffAna/BDTWeights/BPw.root')
shutil.copy2(outFiles[1], '../gd/braa_nohbhe_trk/Bs/EffAna/BDTWeights/Bsw.root')
shutil.copy2(outFiles[1], '../gd/braa_mergebin/Bs/EffAna/BDTWeights/Bsw.root')
