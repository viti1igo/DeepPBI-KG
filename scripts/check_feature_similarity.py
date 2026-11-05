"""Recheck feature_similarity_report.md numbers. Run from the repo root.

Compares the scaler used in code/visualize_features.py (wgs_scaler for both
feature sets) with what DeepPBI-KG.py actually feeds its models (each feature
set scaled by its own scaler).
"""
import pickle,numpy as np,pandas as pd,warnings
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
warnings.filterwarnings('ignore')
for sp,fw,fk,sl in [('Phage','phage_wgs','phage_key_gene',slice(0,1338)),('Bacterium','bacterium_wgs','bacterium_key_gene',slice(1338,None))]:
  w=pd.read_csv(f'template/feature_file/{fw}_feature.csv',index_col=0);k=pd.read_csv(f'template/feature_file/{fk}_feature.csv',index_col=0)
  c=w.index.intersection(k.index);print(sp,w.shape,k.shape,'common',len(c),'same cols',list(w.columns)==list(k.columns))
  w=w.loc[c].values;k=k.loc[c].values
  for name,sw,sk in [('report(wgs scaler both)','wgs','wgs'),('classifier(own scaler)','wgs','key_gene')]:
    S={n:pickle.load(open(f'model/{n}_scaler.pkl','rb')) for n in {sw,sk}}
    a=(w-S[sw].mean_[sl])/S[sw].scale_[sl];b=(k-S[sk].mean_[sl])/S[sk].scale_[sl]
    A=a/np.linalg.norm(a,axis=1,keepdims=True);B=b/np.linalg.norm(b,axis=1,keepdims=True);M=A@B.T
    n=len(c);d=np.diag(M);off=M[~np.eye(n,dtype=bool)]
    p=PCA(2).fit(np.vstack([a,b]));pa,pb=p.transform(a)[:,0],p.transform(b)[:,0]
    print(f' {name}: scaler_len={len(S[sw].mean_)} diag mean={d.mean():.3f} min={d.min():.3f} max={d.max():.3f}; off mean={off.mean():.3f} max={off.max():.3f}; rows own best={(M.argmax(1)==np.arange(n)).sum()}/{n}; cols own best={(M.argmax(0)==np.arange(n)).sum()}/{n}; PC1 r={pearsonr(pa,pb)[0]:.4f} p={pearsonr(pa,pb)[1]:.1e} var={p.explained_variance_ratio_.round(3)}; mean wgs-kg PC1 gap={(pb-pa).mean():.1f}, wgs PC1 sd={pa.std():.1f}')
    print('  mean pearson',np.mean([pearsonr(a[i],b[i])[0] for i in range(n)]).round(3))
