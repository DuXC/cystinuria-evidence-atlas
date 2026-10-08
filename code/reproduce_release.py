"""Recompute atlas v1.0.1 key results using distributed derived TSVs, standard library only.

Usage: python3 code/reproduce_release.py [unpacked_release_directory]
This verifies the distributed data and reports counts; it does not recreate the
original genomic mapping or reassess the underlying experimental images.
"""
import csv,hashlib,json,sys
from pathlib import Path
R=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
S=R/'supplementary_data'
def read(name):
    with (S/name).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
v=read('S1_variant_atlas_488.tsv');f=read('S2_functional_results_105.tsv')
assert len(v)==len({x['variation_id'] for x in v})==488
assert len(f)==len({x['evidence_id'] for x in f})==105
p=[x for x in f if x['counts_as_primary_human_protein']=='True']
t=[x for x in f if x['counts_as_primary_human_transport']=='True']
k=lambda x:(x['gene'],x['protein_variant'])
vk=lambda x:(x['gene'],x['uniprot_change'])
pk={k(x) for x in p};tk={k(x) for x in t}
np=sum(vk(x) in pk for x in v);nt=sum(vk(x) in tk for x in v)
assert (len(p),len(t),len(pk),len(tk),np,nt)==(51,43,25,22,22,20)
assert all((x['has_human_primary_fulltext_function']=='True')==(vk(x) in pk) for x in v)
assert all((x['has_human_primary_fulltext_transport']=='True')==(vk(x) in tk) for x in v)
loo={}
for x in read('S14_source_dependency.tsv'):
    fam=x['omitted_family'];keys={k(z) for z in p if z['experimental_family_id']!=fam}
    retained=sum(vk(z) in keys for z in v);assert retained==int(x['retained_protein_vcv'])
    keys_t={k(z) for z in t if z['experimental_family_id']!=fam}
    assert sum(vk(z) in keys_t for z in v)==int(x['retained_transport_vcv'])
    loo[fam]=retained
dist=['min_arg_heavyatom_A','min_heteromer_heavyatom_A','min_higher_order_heavyatom_A']
thresholds={};members=read('S19_contact_membership.tsv')
for cut in [4,5,6]:
    selected=[x for x in v if x['clinical_group'] in ['VUS','Conflict'] and vk(x) not in pk and any(x[c] and float(x[c])<=cut for c in dist)]
    thresholds[cut]=len(selected)
    assert {x['variation_id'] for x in selected}=={x['variation_id'] for x in members if int(x['contact_cutoff_A'])==cut}
    for row in read('S15_current_contact_thresholds.tsv'):
        if int(row['contact_cutoff_A'])==cut:
            sub=[x for x in selected if row['gene']=='ALL' or x['gene']==row['gene']]
            assert len(sub)==int(row['unmatched_vcv'])
            assert len({(x['gene'],x['uniprot_pos']) for x in sub})==int(row['unique_gene_residues'])
assert list(thresholds.values())==[42,52,58]
byid={x['evidence_id']:x for x in f}
paired=read('S16_paired_assay_summaries.tsv')
assert len(paired)==12
for row in paired:
    for col,value in row.items():assert byid[row['evidence_id']][col]==value,(row['evidence_id'],col)
for row in read('S18_coverage_by_scope.tsv'):
    sub=[x for x in v if row['stratum']=='overall' or (x['gene'] if row['stratum']=='gene' else x['clinical_group'])==row['group']]
    assert len(sub)==int(row['vcv_denominator'])
    assert sum(vk(x) in pk for x in sub)==int(row['protein_vcv'])
    assert sum(vk(x) in tk for x in sub)==int(row['transport_vcv'])
index=read('TABLE_INDEX.tsv');assert len(index)==19
for row in index:
    assert hashlib.sha256((S/row['file']).read_bytes()).hexdigest()==row['output_sha256'],row['file']
result={'status':'PASS','variant_records':len(v),'curated_rows':len(f),'protein_records':np,'transport_records':nt,
 'numeric_or_censored_records':sum(x['has_human_primary_fulltext_summary_numeric']=='True' for x in v),
 'protein_family_substitution_cells':len({(x['experimental_family_id'],*k(x)) for x in p}),
 'transport_family_substitution_cells':len({(x['experimental_family_id'],*k(x)) for x in t}),
 'leave_one_family_out':loo,'current_contact_records':thresholds,'source_summary_rows_verified':len(paired),'table_hashes_verified':len(index)}
assert result['numeric_or_censored_records']==9
print(json.dumps(result,ensure_ascii=False,indent=2))
