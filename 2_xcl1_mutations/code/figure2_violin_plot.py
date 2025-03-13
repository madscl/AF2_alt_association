'''
Creates violin plot for Figure 2
'''

import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd


column_name = ['variant', 'model', 'seed', 'TM-score','reference']
l1=[]
l2=[]
l3=[]
l4=[]
l5=[]
models= 5
key_num = list(range(0, 56, 3))
for m in range(1,models+1): 
    seeds = 5
    for s in range(seeds):
        path1 = './2_xcl1_mutations/TM-scores/recycles_3/only_TMscores_1j9o_model'+str(m)+'_seed'+str(s)+'_3.txt'
        path2 = './2_xcl1_mutations/TM-scores/recycles_3/only_TMscores_2jp1_model'+str(m)+'_seed'+str(s)+'_3.txt'
        with open(path1) as file:
            content = file.readlines()
            # Goes to each line that has a name
            for num in key_num:
                l1.append(content[num].strip())
                l2.append(m)
                l3.append(s)
                l4.append(float(content[num+1].strip()))
                l5.append('chemokine')
        with open(path2) as file:
            content = file.readlines()
            # Goes to each line that has a name
            for num in key_num:
                l1.append(content[num].strip())
                l2.append(m)
                l3.append(s)
                l4.append(float(content[num+1].strip()))
                l5.append('dimer')
df1 = pd.DataFrame(list(zip(l1, l2, l3, l4, l5)), columns=column_name)
print(df1)

plt.figure(figsize=(20,6))

ax=sns.violinplot(data=df1, x="variant",y="TM-score",hue="reference", palette=["red","gold"],scale="width", width=0.5, linewidth=0.5)

for i in range(0,10):
    ax.axvspan(xmin=(-0.5+2*i), xmax=(0.5+2*i),facecolor="lightgrey", alpha=0.3)
ax.set_xlabel('Variant', fontsize=22) 
ax.set_ylabel('TM-score', fontsize=22)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.tick_params(axis='x', labelsize=22) 
ax.tick_params(axis='y', labelsize=22) 
ax.hlines(0.5,xmin=-0.5, xmax=18.5,colors='grey', linestyles='dashed')
ax.margins(x=0)
plt.legend(title='reference',loc='lower left',ncol=4,fontsize=18,title_fontsize=22)
plt.show()

