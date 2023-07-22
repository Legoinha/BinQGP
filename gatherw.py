#! /bin/python
import ROOT as r
import shutil
import re
# ptlist = [5, 7, 10, 15, 20, 50, 60]
# ptlist = [5, 7, 10, 15, 20, 60]

# outFiles = ['BPw.root', 'Bsw.root']

# fout = r.TFile(outFiles[0], 'recreate')
# for pti, ptf in zip(ptlist[:-1], ptlist[1:]):
#     fin = r.TFile(f'./results/Bu/{pti}_{ptf}/mc_validation_plots/weights/weights.root')
#     h = fin.Get(f'weights_BDT_pt_{pti}_{ptf}')
#     h2 = h.Clone()
#     fout.cd()
#     h2.Write()
# fout.Close()

# ptlist = [7, 10, 15, 20, 50]
# fout = r.TFile(outFiles[1], 'recreate')
# for pti, ptf in zip(ptlist[:-1], ptlist[1:]):
#     fin = r.TFile(f'./results/Bs/{pti}_{ptf}/mc_validation_plots/weights/weights.root')
#     h = fin.Get(f'weights_BDT_pt_{pti}_{ptf}')
#     h2 = h.Clone()
#     fout.cd()
#     h2.Write()
# fout.Close()

# shutil.copy2(outFiles[0], '../gd/braa_nohbhe_trk/BP/EffAna/BDTWeights/BPw.root')
# shutil.copy2(outFiles[0], '../gd/braa_mergebin_git/BP/EffAna/BDTWeights/BPw.root')
# shutil.copy2(outFiles[1], '../gd/braa_nohbhe_trk/Bs/EffAna/BDTWeights/Bsw.root')
# shutil.copy2(outFiles[1], '../gd/braa_mergebin_git/Bs/EffAna/BDTWeights/Bsw.root')


bpptlist = [5, 15, 60]
bsptlist = [7, 15, 50]
combptlist = (bpptlist, bsptlist)
varlist = ['BsvpvSig', 'Bchi2cl', 'Bmumumass']
btypes = ['BP', 'Bs']
btypes_low = ['Bp', 'Bs']
cats = ['Bu', 'Bs']

# for bt, bcat, ptlist in zip(btypes, cats, combptlist):
#     for pti, ptf in zip(ptlist[:-1], ptlist[1:]):
#         for v in varlist:
#             # weight_file = f"./results/{bcat}/inclusive_{pti}_{ptf}_3gauss/mc_validation_plots/weights/weights.root"
#             # shutil.copy2(weight_file, f'../gd/braa_nohbhe_trk/{bt}/EffAna/BDTWeights/{bt}weights.root')
#             fcopy = f'./results/{bcat}/inclusive_{pti}_{ptf}_3gauss/mc_validation_plots/mc_sp/pdfs/{v}_mc_validation_{bcat}.pdf'
#             shutil.copy2(fcopy, f'../overleaf/Plots/Syst/MCData/{bt}/{v}_pt_inclusive_{pti}_{ptf}.pdf')
#                 # ystr = "_y15-24" if pti < 10 else ""
#                 # fcopy = f'./results/{bcat}/{pti}_{ptf}/mc_validation_plots/mc_sp/pdfs/{v}_mc_validation_{bcat}{ystr}.pdf'
#                 # print(fcopy)

#                 # shutil.copy2(fcopy, f'../overleaf/Plots/Syst/MCData/BP/{v}_pt_{pti}_{ptf}.pdf')

#         # fcopy = f'/home/tasheng/braa/BinQGP_trk5/results/Bu/{pti}_{ptf}/mc_validation_plots/mc_sp/pdfs/BDT_pt_{pti}_{ptf}_mc_validation_Bu{ystr}.pdf'
#         # print(fcopy)
#         # shutil.copy2(fcopy, f'../overleaf/Plots/Syst/MCData/BP/')


varlist = ['BsvpvDistRatio', 'Bchi2cl']
outFiles = ['BPweights.root', 'Bsweights.root']
suffix = ['Low', 'High']
for bt, bcat, ptlist, outstr in zip(btypes, cats, combptlist, outFiles):
    for pti, ptf, suf in zip(ptlist[:-1], ptlist[1:], suffix):
        outfile = outstr[:-5] + suf + '.root'
        fout = r.TFile(outfile, 'recreate')
        rep = '_rep' if bt == 'BP' else ''
        for v in varlist:
            infile = f'./results/{bcat}/inclusive_{pti}_{ptf}_{v}{rep}/mc_validation_plots/weights/weights.root'
            # infile = f'./results/{bcat}/inclusive_{pti}_{ptf}/mc_validation_plots/weights/weights.root'

            fin = r.TFile(infile)
            # print(infile)
            if v == 'BsvpvDistRatio':
                oth = 'Bchi2cl'
            elif v == 'Bchi2cl':
                oth = 'BsvpvSig'
            else:
                print('error')
            for cate in ['mc', 'data']:
                h = fin.Get(f'{cate}_{oth}')
                h2 = h.Clone(f'{cate}_{oth}')
                # print(h.GetBinContent(1))
                fout.cd()
                h2.Write()

            fcopy = f'./results/{bcat}/inclusive_{pti}_{ptf}_{v}{rep}/mc_validation_plots/mc_sp/pdfs/{oth}_mc_validation_{bcat}.pdf'
            shutil.copy2(fcopy, f'../overleaf/Plots/Syst/MCData/{bt}/{oth}_pt_inclusive_{pti}_{ptf}.pdf')
        fout.Close()
        shutil.copy2(outfile, f'../gd/braa_nohbhe_trk/{bt}/EffAna/BDTWeights/')
        shutil.copy2(outfile, f'../gd/braa_mergebin_git/{bt}/EffAna/BDTWeights/')
