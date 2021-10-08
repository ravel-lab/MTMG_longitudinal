import pandas as pd

MT_cst_I_IV = pd.read_csv("CST_I_IV_kegg_module.csv",sep=",",index_col=0)
MT_cst_III_IV = pd.read_csv("CST_III_IV_kegg_module.csv",sep=",",index_col=0)
MT_cst_I_III = pd.read_csv("CST_I_III_kegg_module.csv",sep=",",index_col=0)

def f(row):
    
    row = row.abs()
    x = row.max()
    
    return x


count =0
for compar in ['I_III','I_IV','III_IV']:
    print(compar)

    if compar == 'I_III':

        df1 = MT_cst_I_III

    elif compar == 'I_IV':

        df1 = MT_cst_I_IV

    elif compar == 'III_IV':

        df1 = MT_cst_III_IV

    df1_fc = df1.drop(['AxeExpr','t','P','adjP','Z','Count'],axis=1)
    df1_count = df1.drop(['logFC','AxeExpr','t','P','adjP','Z'],axis=1)

    df1_fc = pd.pivot_table(data=df1_fc,index="Pathway",columns="Sign",values='logFC')
    df1_count = pd.pivot_table(data=df1_count,index="Pathway",columns="Sign",values='Count')

    df1_fc = df1_fc.fillna(0)
    df1_count = df1_count.fillna(0)

    df1_fc['maxFC'] = df1_fc.apply(f,axis=1)
    df1_count['maxCount'] = df1_count.apply(f,axis=1)
    
    df1_fc_trim = df1_fc[df1_count['maxCount']>=1]
    df1_fc_trim = df1_fc_trim[df1_fc_trim['maxFC'] >= 4]

    df1_fc_trim = df1_fc_trim.drop(['maxFC'],axis=1)

    df1_fc_trim.columns = [compar.split("_")[1],compar.split("_")[0]]

    if count==0:

        df_out_fc = df1_fc_trim
        count+=1

    elif count >0:
        
        df_out_fc = pd.merge(left=df_out_fc,right=df1_fc_trim,left_index=True,right_index=True,how="outer")

df_out_fc.to_csv("MT_module_summary_aveFC4.csv",sep=",")



