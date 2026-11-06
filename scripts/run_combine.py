import os,sys
import numpy as np
from optparse import OptionParser
def readNorm(f_cscCard):
    f = open(f_cscCard,"r")
    # norm = float(f.readline().split()[3])
    Lines = f.readlines()
    for line in Lines:
        l  = line.split()
        if len(l) == 0:continue
        if l[0] == 'rate':
            return float(l[1])
def insert(originalfile,string):
    with open(originalfile,'r') as f:
        with open('newfile.txt','w') as f2:
            f2.write(string)
            f2.write(f.read())
    os.rename('newfile.txt',originalfile)

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--dryRun', dest='dryRun', action="store_true",default = False, help='print cmds')
    (options, args) = parser.parse_args()

    #list all the DT cards
    CMSSW_BASE = os.getenv('CMSSW_BASE')

    # ---- FILLED IN ---------------------------------------------------------
    # create_datacard.ipynb writes the cards with
    #     outDataCardsDir = os.getcwd().replace('scripts', 'combine/datacards/')
    # i.e. $CMSSW_BASE/src/MDS_CMSDAS/combine/datacards/ (one .txt per
    # ggH_HToSSTobbbb_mh125_mx40_ctau<ct>mm point).
    # os.getenv('CMSSW_BASE') has NO trailing slash, so the string appended here
    # must start with '/'. It must also keep the directory name 'datacards',
    # because output_tree_dir below is built with .replace('datacards','limitTrees').
    card_dir= CMSSW_BASE + "/src/MDS_CMSDAS/combine/datacards/"

    output_tree_dir = card_dir.replace('datacards', 'limitTrees')
    print(output_tree_dir)
    if not os.path.exists(output_tree_dir):os.system("mkdir -p %s"%output_tree_dir)
    ctau_list = [5,10,15,20,30,40, 50,60, 100, 125,150,200,300,500,600,700,800,900,1000, 2000,3000,4000, 5000, 6000,7000,8000,10000, 20000,30000,50000, 100000, 200000, 300000, 500000, 1000000, 2000000, 3000000, 5000000, 6000000, 10000000] #mm



    ##ctau_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    # 'sample_name' must match the modelName built in create_datacard.ipynb:
    # 'ggH_HToSSTo' + decay + '_mh125_mx' + str(m) + '_ctau' + ctau + 'mm'.
    # For the Exercise (mS = 15 and 55 GeV) change mx40 -> mx15 / mx55 here and
    # re-run, after re-running create_datacard.ipynb with those masses added.
    sample_name = 'ggH_HToSSTobbbb_mh125_mx40'
    for ctau in ctau_list:
        ctau_string = str(ctau)
        name = "{sample}_ctau{ctau}mm".format(sample=sample_name,ctau=ctau_string)

        ### check all files exist
        if not os.path.isfile(card_dir + name + '.txt'):
    	    print(card_dir + name + '.txt')
    	    continue
        #normalization = readNorm(twoTag_card_dir + twoTag_name + '.txt')
        norm = readNorm(card_dir + name + '.txt')

        if norm == 0.0:continue
        # ---- FILLED IN -------------------------------------------------
        # readNorm() returns the SECOND token of the 'rate' line of the card, which
        # helper.make_datacard() writes as  rate  <signal yield of chA> 1 <sigB> 1 ...
        # so 'norm' here is the signal yield in the signal region (bin A), sigA.
        # The card is scaled with the 'norm' rateParam, so we run Combine with
        #   norm = 1./sigA
        # (the rule of thumb quoted in the episode). The fitted signal strength r is
        # then O(1-15) for every lifetime; the scale factor is multiplied back when the
        # limits are plotted (limitPlot.ipynb divides the limit by sigA).
        #
        # The original string had two problems: '--freezeParameter=norm' is not a
        # Combine option (it is '--freezeParameters <name>', plural, space separated),
        # and '{card_dir}/{card}.txt' doubled the slash now that card_dir ends in '/'.
        comb_limit = "combine -M AsymptoticLimits {card_dir}{card}.txt -n.{card} --setParameters norm={norm} --freezeParameters norm".format(card_dir=card_dir,card=name,norm = 1./norm)
        #print(comb_limit)
        os.system(comb_limit)
        os.system("mv higgsCombine.{}.AsymptoticLimits.mH120.root {}/".format(name, output_tree_dir))
