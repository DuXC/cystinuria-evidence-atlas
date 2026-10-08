"""Reproducible descriptive figures for CBC manuscript v7."""
from pathlib import Path
import hashlib,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RELEASE_MODE=(ROOT/'supplementary_data').is_dir()
DATA=ROOT/'analysis_source_data' if RELEASE_MODE else ROOT/'data/processed/editorial_20261008'
OUT=ROOT/'figures' if RELEASE_MODE else ROOT/'deliverables/editorial_revision_20261008_v6_1/figures'
OLD=None if RELEASE_MODE else ROOT/'deliverables/author_review_20261008_v5_5/figures'
LOG=ROOT/'reproduction' if RELEASE_MODE else ROOT/'logs/editorial_20261008_v6_1'
LOG.mkdir(parents=True,exist_ok=True)
C={'navy':'#17324D','teal':'#007C83','orange':'#C77821','muted':'#516273','pale':'#E6ECEF'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':7.5,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.6,'axes.edgecolor':'#87939E','text.color':C['navy'],'axes.labelcolor':C['navy'],
 'xtick.color':C['muted'],'ytick.color':C['navy'],'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,
 'savefig.facecolor':'white'})
OUT.mkdir(parents=True,exist_ok=True)
provenance=[]
def save(fig,name,source):
    for ext in ['png','svg','pdf']:
        p=OUT/f'{name}.{ext}';fig.savefig(p,dpi=600)
        provenance.append({'file':p.name,'source':source,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'mode':'NEW_DESCRIPTIVE_PLOT'})
    fig.savefig(LOG/f'{name}_preview.png',dpi=150)
    plt.close(fig)
def label(ax,letter,title):
    ax.text(-.13,1.10,letter,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom')
    ax.set_title(title,loc='left',pad=17,fontweight='bold')
def clear_y(ax):
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)

def coverage():
    fig=plt.figure(figsize=(183/25.4,155/25.4))
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.1],left=.21,right=.96,top=.88,bottom=.09,hspace=.72,wspace=1.0)
    a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,:])
    vals=np.array([22,20,9,1]);yy=np.arange(4)
    a.barh(yy,vals/488*100,color=[C['teal']]*3+[C['orange']],height=.46)
    for y,n in zip(yy,vals):a.text(n/488*100+.13,y,f'{n}/488',va='center',fontsize=7.5)
    a.set_yticks(yy,['Protein','Transport','Numeric / bound','Allele-specific RNA'])
    a.invert_yaxis();a.set_xlim(0,6.3);a.set_xticks([0,2,4,6]);a.set_xlabel('Records with evidence (%)');clear_y(a)
    label(a,'a','Coverage by endpoint')
    cl=pd.read_csv(DATA/'coverage_by_scope.tsv',sep='\t').query("stratum == 'classification'")
    b.scatter(cl.protein_percent,np.arange(5),s=23,color=C['navy'])
    for y,(_,row) in enumerate(cl.iterrows()): b.text(row.protein_percent+2,y,f'{row.protein_vcv}/{row.vcv_denominator}',va='center',fontsize=7.5)
    b.set_yticks(np.arange(5),['P/LP','VUS','Conflicting','B/LB','Other']);b.invert_yaxis();b.set_xlim(0,56);b.set_xticks([0,20,40]);b.set_xlabel('Protein-matched records (%)');clear_y(b)
    label(b,'b','Coverage by archive label')
    loo=pd.read_csv(DATA/'leave_one_family_out.tsv',sep='\t');yy=np.arange(6)
    for y,n in zip(yy,loo.retained_protein_vcv):c.plot([n,22],[y,y],color=C['pale'],lw=4,zorder=1)
    c.axvline(22,color=C['muted'],lw=.8,ls='--',zorder=0)
    c.scatter(loo.retained_protein_vcv,yy,s=32,color=[C['orange'] if n<22 else C['teal'] for n in loo.retained_protein_vcv],zorder=3)
    for y,(_,row) in enumerate(loo.iterrows()):c.text(22.7,y,f'{row.retained_protein_vcv} retained   ({row.lost_protein_vcv} lost)',va='center',fontsize=8)
    c.set_yticks(yy,['Font 2001','Pineda 2004','Shigeta 2006','Bartoccioni 2008','Yan 2020','He 2026']);c.invert_yaxis();c.set_xlim(13.5,29);c.set_xticks([14,16,18,20,22]);c.set_xlabel('Protein-matched VCV records after removing one family');clear_y(c)
    label(c,'c','Dependence on the checked source families')
    fig.text(.045,.965,'Functional coverage is limited and source-dependent',fontsize=12,fontweight='bold',va='top')
    save(fig,'Figure1_coverage_and_source_dependence','S14_source_dependency.tsv; S18_coverage_by_scope.tsv; S1_variant_atlas_488.tsv')

def assay():
    d=pd.read_csv(DATA/'paired_assay_summaries.tsv',sep='\t')
    fig=plt.figure(figsize=(183/25.4,150/25.4))
    gs=fig.add_gridspec(1,2,width_ratios=[1.3,1],left=.12,right=.96,top=.79,bottom=.29,wspace=.30)
    a=fig.add_subplot(gs[0]);b=fig.add_subplot(gs[1]);b.axis('off')
    variants=['A70V','A182T','G105R','R333W','V170M','A354T']
    for y,v in enumerate(variants):
        rows=d[d.protein_variant==v];f=rows[rows.study_id=='Font2001'].iloc[0];h=rows[rows.study_id=='He2026'].iloc[0]
        a.scatter(f.mean_percent_WT,y-.12,s=27,color=C['navy'],label='Font 2001' if y==0 else None,zorder=4)
        if h.numeric_kind=='LEFT_CENSORED':
            a.plot([0,h.upper_percent_WT],[y+.14]*2,color=C['teal'],lw=1.3)
            a.scatter(h.upper_percent_WT,y+.14,marker='<',facecolors='white',edgecolors=C['teal'],s=35,zorder=4)
            a.text(8,y+.17,'<5',fontsize=7.5,color=C['teal'],va='center')
        else:
            a.errorbar(h.mean_percent_WT,y+.14,xerr=h.sd_percent_WT,fmt='s',color=C['teal'],ms=4,capsize=2,lw=1.1,label='He 2026: mean ± SD' if y==0 else None)
    a.axvline(100,color=C['pale'],ls='--',lw=1)
    a.set_yticks(np.arange(6),variants);a.set_ylim(5.55,-.65);a.set_xlim(-3,107);a.set_xticks([0,25,50,75,100]);a.set_xlabel('Activity (% of each study’s wild type)');clear_y(a)
    a.legend(frameon=False,fontsize=7.8,loc='lower left',bbox_to_anchor=(-.1,1.02),handlelength=1.5)
    a.text(-.22,1.21,'a',transform=a.transAxes,fontsize=11,fontweight='bold');a.text(0,1.21,'Same substitution, different assays',transform=a.transAxes,fontsize=9,fontweight='bold')
    b.text(0,1.21,'b',transform=b.transAxes,fontsize=11,fontweight='bold');b.text(.1,1.21,'R365W: culture condition',transform=b.transAxes,fontsize=9,fontweight='bold')
    b.text(0,1.04,'Pineda 2004 · SLC3A1',fontsize=8.5,fontweight='bold')
    b.text(0,.84,'33 °C culture',fontsize=10,fontweight='bold',color=C['teal'])
    b.text(0,.64,'49 ± 4% WT',fontsize=17,fontweight='bold',color=C['teal'])
    b.text(0,.51,'Mean ± SEM\n11 independent experiments',fontsize=8.5,linespacing=1.5)
    b.plot([0,1],[.42,.42],transform=b.transAxes,color=C['pale'],lw=1)
    b.text(0,.30,'37 °C culture',fontsize=10,fontweight='bold')
    b.text(0,.13,'No uptake above background',fontsize=8.3,fontweight='bold')
    b.text(0,-.02,'Qualitative source observation',fontsize=8,color=C['muted'])
    b.set_xlim(0,1);b.set_ylim(0,1)
    fig.text(.035,.965,'Functional readouts retain their experimental context',fontsize=12,fontweight='bold',va='top')
    fig.text(.12,.135,'Font: 20 µM radiolabeled cystine · HeLa cells\nHe: 200 µM selenocystine · HEK293 cells · 30 min\nHe: 3 independent transfections, technical triplicates',fontsize=8,linespacing=1.5)
    fig.text(.625,.135,'R365W uptake measured at 37 °C.\nCulture and measurement temperatures\nare distinct experimental conditions.',fontsize=8,linespacing=1.5)
    fig.text(.12,.055,'Source summaries are displayed separately; differences do not isolate a system, substrate or timing effect.',fontsize=7.5,color=C['muted'])
    save(fig,'Figure3_assay_context','S16_paired_assay_summaries.tsv; S2_functional_results_105.tsv R2P001/R2P002')

def thresholds():
    d=pd.read_csv(DATA/'current_contact_thresholds.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(140/25.4,90/25.4));fig.subplots_adjust(left=.16,right=.96,top=.78,bottom=.18)
    x=np.arange(3);bottom=np.zeros(3)
    for gene,color in [('SLC3A1',C['navy']),('SLC7A9',C['teal'])]:
        vals=d[d.gene==gene].unmatched_vcv.to_numpy();ax.bar(x,vals,bottom=bottom,width=.5,color=color,label=gene)
        for xx,v,bt in zip(x,vals,bottom):ax.text(xx,bt+v/2,str(v),ha='center',va='center',color='white',fontsize=10)
        bottom+=vals
    for xx,n in zip(x,bottom):ax.text(xx,n+1.4,str(int(n)),ha='center',fontweight='bold',fontsize=11)
    ax.set_xticks(x,['4 Å','5 Å','6 Å']);ax.set_ylim(0,68);ax.set_ylabel('VCV records without a checked protein result');ax.set_xlabel('Contact distance threshold');ax.legend(frameon=False,ncol=2,loc='lower left',bbox_to_anchor=(-.05,1.01))
    fig.suptitle('Current contact-region candidates depend on the cutoff',fontsize=11,fontweight='bold',x=.07,ha='left',y=.98)
    save(fig,'FigureS3_current_contact_thresholds','S15_current_contact_thresholds.tsv; S19_contact_membership.tsv')

def main():
    coverage();assay();thresholds()
    if RELEASE_MODE:
        (LOG/'figure_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
        print('Recreated Figure 1, Figure 3 and Supplementary Figure S3.')
        return
    for src,dst in [('Figure2_experimental_evidence_matrix','Figure2_experimental_evidence_matrix'),('Figure1_variant_sequence_atlas','FigureS1_variant_sequence_atlas'),('Figure3_original_contact_region_cohort','FigureS2_original_contact_region_cohort')]:
        for ext in ['png','svg','pdf']:
            a=OLD/f'{src}.{ext}';b=OUT/f'{dst}.{ext}';shutil.copy2(a,b)
            provenance.append({'file':b.name,'source':str(a.relative_to(ROOT)),'sha256':hashlib.sha256(b.read_bytes()).hexdigest(),'mode':'BYTE_IDENTICAL_REUSE'})
    (LOG/'figure_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Created 3 new descriptive figures and retained 3 validated figures; PNG/SVG/PDF each.')
if __name__=='__main__': main()
