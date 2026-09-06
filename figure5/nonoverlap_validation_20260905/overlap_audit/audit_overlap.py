from pathlib import Path
import csv, gzip, re, json, hashlib
import numpy as np

OUT=Path(__file__).resolve().parent
BASE=Path('G:/1Yunvjian/periodontitis基因/GEO')
ARCH=Path('G:/1Yunvjian/1A玉女煎/code/coredata/bulk')
def read_matrix(path):
    meta={}; values=[]; probes=[]
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            row=next(csv.reader([line],delimiter='\t'))
            if not row: continue
            if row[0]=='!series_matrix_table_begin': break
            if row[0].startswith('!Sample_'): meta.setdefault(row[0],[]).append(row[1:])
        header=next(csv.reader([next(f)],delimiter='\t'))
        for line in f:
            if line.startswith('!series_matrix_table_end'): break
            row=next(csv.reader([line],delimiter='\t'))
            probes.append(row[0]); values.append([float(v) for v in row[1:]])
    records=[]
    for gsm,title in zip(meta['!Sample_geo_accession'][0],meta['!Sample_title'][0]):
        m=re.search(r'sample (\d+) from patient (\d+), (.*)',title,re.I)
        assert m,title
        site,patient,phenotype=m.groups()
        records.append(dict(sample_id=gsm,patient=int(patient),site=int(site),group='NC' if 'unaffected' in phenotype.lower() else 'PD',title=title))
    return records,probes,np.array(values,dtype=float)
def write_csv(name, rows):
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

cache=OUT/'matrix_cache.npz'
if cache.exists():
    z=np.load(cache,allow_pickle=True);a=z['a'].tolist();b=z['b'].tolist();pa=z['pa'].tolist();pb=pa;xa=z['xa'];xb=z['xb']
else:
    a,pa,xa=read_matrix(BASE/'数据处理/GSE10334_series_matrix.txt.gz')
    b,pb,xb=read_matrix(BASE/'外部数据/GSE16134_analysis/GSE16134_series_matrix.txt.gz')
    np.savez(cache,a=a,b=b,pa=pa,xa=xa,xb=xb)
assert pa==pb
split=list(csv.DictReader((ARCH/'patient_split.csv').open(encoding='utf-8-sig')))
dev={r['sample_id']:r for r in split if r['dataset']=='GSE10334'}
bykey={(r['patient'],r['site']):(j,r) for j,r in enumerate(b)}
assert len(bykey)==len(b)
# All-probe nearest-neighbour correlation is independent of the title match.
za=xa-xa.mean(axis=0);zb=xb-xb.mean(axis=0)
za/=np.linalg.norm(za,axis=0);zb/=np.linalg.norm(zb,axis=0)
corr=za.T@zb
rows=[]
for i,r in enumerate(a):
    ranked=np.argsort(corr[i]);nearest=int(ranked[-1])
    j=nearest; s=b[j]
    if (r['patient'],r['site']) != (s['patient'],s['site']):
        print('TITLE DISCREPANCY',r,s,'correlation',float(corr[i,j]),flush=True)
    rows.append(dict(GSE10334_sample=r['sample_id'],GSE16134_sample=s['sample_id'],patient=r['patient'],site=r['site'],group=r['group'],external_group=s['group'],used_in_current_development=r['sample_id'] in dev,development_partition=dev.get(r['sample_id'],{}).get('cohort',''),development_patient_key=dev.get(r['sample_id'],{}).get('patient_id',''),external_patient=s['patient'],external_site=s['site'],title_key_match=(r['patient'],r['site'])==(s['patient'],s['site']),all_probe_correlation=float(corr[i,j]),nearest_expression_sample=b[nearest]['sample_id'],nearest_matches_title=(r['patient'],r['site'])==(s['patient'],s['site']),second_best_correlation=float(corr[i,ranked[-2]]),max_absolute_expression_difference=float(np.max(np.abs(xa[:,i]-xb[:,j])))))
devpatients={r['patient'] for r in rows if r['used_in_current_development']}
originalpatients={r['patient'] for r in a}
dev_external={r['GSE16134_sample'] for r in rows if r['used_in_current_development'] and r['title_key_match']}
disposition=[]
for r in b:
    disposition.append(dict(**r,exact_sample_in_current_development=r['sample_id'] in dev_external,patient_in_current_development=r['patient'] in devpatients,patient_in_original_GSE10334=r['patient'] in originalpatients,recommended_action='exclude_shared_patient' if r['patient'] in originalpatients else 'retain_for_reassessment'))
retained=[r for r in disposition if r['recommended_action']=='retain_for_reassessment']
summary=dict(original_GSE10334_samples=len(a),original_GSE10334_patients=len(originalpatients),GSE16134_samples=len(b),GSE16134_patients=len({r['patient'] for r in b}),common_probes=len(pa),current_GSE10334_development_samples=len(dev),current_GSE10334_development_patients=len(devpatients),development_training=sum(r['development_partition']=='Training' for r in rows),development_test=sum(r['development_partition']=='Internal test' for r in rows),all_original_samples_title_mapped=sum(r['title_key_match'] for r in rows),all_original_samples_expression_nearest_match=sum(r['nearest_matches_title'] for r in rows),minimum_matched_correlation=min(r['all_probe_correlation'] for r in rows if r['title_key_match']),maximum_matched_absolute_difference=max(r['max_absolute_expression_difference'] for r in rows if r['title_key_match']),external_samples_from_development_patients=sum(r['patient_in_current_development'] for r in disposition),retained_samples=len(retained),retained_patients=len({r['patient'] for r in retained}),retained_PD=sum(r['group']=='PD' for r in retained),retained_NC=sum(r['group']=='NC' for r in retained),retained_patient_ids=sorted({r['patient'] for r in retained}))
assert len(dev)==sum(r['used_in_current_development'] for r in rows)
# Unresolved matches are exported for review.
# Patient discrepancies require further review.
# Phenotype discrepancies require further review.
summary.update(current_development_samples_with_supported_external_match=len(dev_external),unresolved_original_samples=[r['GSE10334_sample'] for r in rows if not r['title_key_match']],unresolved_current_development_samples=[r['GSE10334_sample'] for r in rows if not r['title_key_match'] and r['used_in_current_development']],retained_NC_patients=len({r['patient'] for r in retained if r['group']=='NC'}),retained_PD_patients=len({r['patient'] for r in retained if r['group']=='PD'}))
for r in rows:
    r['match_evidence']='matching_patient_site_phenotype_and_all_probe_nearest_neighbour' if r['title_key_match'] and r['group']==r['external_group'] else 'unresolved_no_identity_claim'
    if not r['title_key_match']: r['GSE16134_sample']=''
write_csv('GSE10334_GSE16134_crosswalk.csv',rows)
write_csv('GSE16134_sample_disposition.csv',disposition)
write_csv('GSE16134_nonoverlapping_patient_samples.csv',retained)
(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
