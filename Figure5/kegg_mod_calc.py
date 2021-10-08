import pandas as pd
data = pd.read_csv("MT_cstIII_IV_kegg_de.csv",sep=",")
data.columns = ['KEGG','logFC','AxeExpr','t','P','adjP','Z']
mod_key = pd.read_csv("~/bioinformatics_tools/kegg/ko/KEGG_MOD_map.txt",sep="\t")
mod_out = pd.merge(left=mod_key,right=data,on="KEGG",how="left")
mod_out = pd.merge(left=mod_key,right=data,on="KEGG",how="right")
mod_group = mod_out.groupby("Module").mean()
mod_out_trim = mod_out[mod_out['adjP'] <0.05]
mod_out_trim = mod_out[mod_out['logFC'].abs() >= 2]
def f(row):
    if row['logFC'] > 0:
        val='Pos'
    if row['logFC'] < 0:
        val='Neg'
    return val

mod_out_trim['Sign'] = mod_out_trim.apply(f,axis=1)
mod_group = mod_out_trim.groupby(["Module",'Sign']).mean()
mod_group_test = mod_out_trim.groupby(["Module",'Sign']).count()
mod_group_test = mod_group_test.drop(['logFC','AxeExpr','t','P','adjP','Z'],axis=1)

mod_group_test=mod_group_test.reset_index()
mod_group = mod_group.reset_index()
mod_group['Count'] = mod_group_test['KEGG']
print(mod_group.head())
mod_group.to_csv("CST_III_IV_kegg_module.csv",sep=",",index=None)