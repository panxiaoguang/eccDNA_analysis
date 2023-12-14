import json
from jinja2 import Environment, PackageLoader
from statistics import median
import plotly
import plotly.figure_factory as ff


env = Environment(loader=PackageLoader('report', 'templates'))
template = env.get_template('index.html')

### get fastp result
def get_fastp(fp):
    fastp_data = json.load(open(fp))
    SeqInfos = []
    names = ['Total Reads','Total Bases','Q20 Bases','Q30 Bases','Q20 Rate','Q30 Rate','Read1 Mean Length','Read2 Mean Length','GC Content']
    SeqInfos.append({"item":"Sequencing Mode","value":fastp_data['summary']['sequencing']})
    for k,v in zip(names,fastp_data['summary']["after_filtering"].values()):
        if isinstance(v,int):
            SeqInfos.append({"item":k,"value":'{:,}'.format(v)})
        else:
            SeqInfos.append({"item":k,"value":round(v,2)})
    return SeqInfos

SeqInfos = get_fastp(snakemake.input['fp'])

### get base stats

def get_ecc_numbers(num):
    with open(num,'r') as f:
        rawEccNumber = int(f.read().strip())
    return rawEccNumber

rawECCNumber = get_ecc_numbers(snakemake.input['eccN'])
rawECCNumber = '{:,}'.format(rawECCNumber)
def get_ecc_length(length):
    lgth = []
    with open(length,'r') as f:
        for line in f:
            lgth.append(int(line.strip()))
    return lgth

lengths = get_ecc_length(snakemake.input['lgth'])
length = median(lengths)
length = '{:,}'.format(length)
def get_gc(gc):
    gcData = []
    with open(gc,'r') as f:
        for line in f:
            gcData.append(float(line.strip()))
    return gcData

gcs = get_gc(snakemake.input['gc'])
GCs = round(median(gcs),2)
### get mapping stats
def get_total_map(mp):
    with open(mp,"r") as f:
        sample,chrom,mt,unmap = f.read().strip().split('\t')
    return sample,int(chrom),int(mt),int(unmap)

sm,chrom,mt,unmap = get_total_map(snakemake.input['mpTotal'])
mapMt = round(mt/(chrom+mt),2)
chrom = '{:,}'.format(chrom)
unmap = '{:,}'.format(unmap)
Mappings = [{"item":"Total Mapped","value":chrom},{"item":"Mapped Ratio of chrM","value":mapMt},{"item":"Unmapped","value":unmap}]
    

### plot ecc length 

def plot_length(length):
    ## only express length <2000
    hist_data = [lg for lg in lengths if lg <2000]
    group_labels = ['eccDNA length']
    fig = ff.create_distplot([hist_data], group_labels, bin_size=0.2,show_rug=False,show_hist=False)
    fig.update_layout(template="simple_white",xaxis_title="length",yaxis_title="density",title="length distribution",height=300)
    fig.update_yaxes(showgrid=True)
    return plotly.io.to_html(fig,full_html=False,include_plotlyjs=True,auto_play=False)

graph = plot_length(lengths)

def plot_GC(gc):
    ## only express length <2000
    hist_data = gc
    group_labels = ['eccDNA GC']
    fig = ff.create_distplot([hist_data], group_labels, bin_size=0.2,show_rug=False,show_hist=False)
    fig.update_layout(template="simple_white",xaxis_title="GC",yaxis_title="density",title="GC distribution",height=300)
    fig.update_yaxes(showgrid=True)
    return plotly.io.to_html(fig,full_html=False,include_plotlyjs=True,auto_play=False)

graphGC = plot_GC(gcs)
### get split and discordant and circle scores

def get_ecc_statistic(ecc):
    split = []
    discordant = []
    circle = []
    with open(ecc,'r') as f:
        for line in f:
            if line.startswith("chrom"):
                continue
            else:
                _,_,_,dc,sp,score,*_ = line.strip().split('\t')
                split.append(int(sp))
                discordant.append(int(dc))
                circle.append(float(score))
    return median(split),median(discordant),median(circle)

split,discordant,circle = get_ecc_statistic(snakemake.input['ecc'])

eccDNAs = [{"item":"Median Split Reads","value":split},{"item":"Median Discordant Reads","value":discordant},{"item":"Median Circle Score","value":circle}]

### get sample stats

Samples = [{"item":"Sample Name","value":snakemake.wildcards.sample},{"item":"Sequencing Platform","value":"CircleSeq"},{"item":"Pipeline Version","value":"0.0.1"},{"item":"Author","value":"Pan Xiaoguang"}]

with open(snakemake.output[0],"w") as f:
    f.write(template.render(rawEccNumber=rawECCNumber,length=length,GCs=GCs,SeqInfos=SeqInfos,Mappings=Mappings,eccDNAs=eccDNAs,Samples=Samples,graph=graph,graphGC=graphGC))
