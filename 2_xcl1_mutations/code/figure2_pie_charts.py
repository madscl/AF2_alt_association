'''
This determines if the TM-score for a given prediction is more chemokine or alternate (if 1 is above 0.5 for one).
If it's below for both then it's considered unstructured
For Figure 2
'''


import matplotlib.pyplot as plt 
import numpy as np 



dimer_file_names = []
chemokine_file_names = []

anc0 = []
anc2 = []
anc2a = []
anc2b = []
anc2c = []
anc2d = []
anc2e = []
anc2f = []
anc2g = []
anc2h = []
anc2i = []
anc2j = []
anc2k = []
anc2l = []
anc2m = []
anc2n = []
anc3 = []
anc4= []
xcl1 = []

models= 5
key_num = list(range(0, 56, 3))
for m in range(1,models+1): 
    # list of variatn names
    l_n = []
    # list of TM scores compared to chemokine
    l_c = []
    # list of TM scores compared to alternate
    l_a = []
    seeds = 5
    for s in range(seeds):
        path1 = './2_xcl1_mutations/TM-scores/recycles_3/only_TMscores_1j9o_model'+str(m)+'_seed'+str(s)+'_3.txt'
        path2 = './2_xcl1_mutations/TM-scores/recycles_3/only_TMscores_2jp1_model'+str(m)+'_seed'+str(s)+'_3.txt'
        with open(path1) as file:
            content = file.readlines()
            # Goes to each line that has a name
            for num in key_num:
                l_n.append(content[num].strip())
                l_c.append(float(content[num+1].strip()))
        with open(path2) as file:
            content = file.readlines()
            # Goes to each line that has a name
            for num in key_num:
                l_n.append(content[num].strip())
                l_a.append(float(content[num+1].strip()))
    for variant in range(0,19):
        # list of more alternate/more chemokine for each protein variant
        l_v = []
        # list for pie chart ratios 
        l_pc = []
        r = 0
        y = 0
        b = 0
        for seed in range(variant, 95, 19):
            a = l_a[seed]
            c = l_c[seed]
            if a > c and a >= 0.5:
                l_v.append('a')
                dimer_file_names.append(l_n[variant]+'_scores_rank_*_alphafold2_model_'+str(m)+'_seed_00'+str(s))
            elif c > a and c >= 0.5:
                l_v.append('c')
                chemokine_file_names.append(l_n[variant]+'_scores_rank_*_alphafold2_model_'+str(m)+'_seed_00'+str(s))
            elif a < 0.5 and c < 0.5:
                l_v.append('n')
        print(l_n[variant], l_v)
        # Determining pie chart ratios
        r = l_v.count('c')
        r = float(r/5)
        y = l_v.count('a')
        y = float(y/5)
        b = l_v.count('n')
        b = float(b/5)
        l_pc = [r, y, b]
        print(l_pc)
        if l_n[variant] == 'anc0':
            anc0.append(l_pc)
        elif l_n[variant] == 'anc2':
            anc2.append(l_pc)
        elif l_n[variant] == 'anc2a':
            anc2a.append(l_pc)
        elif l_n[variant] == 'anc2b':
            anc2b.append(l_pc)
        elif l_n[variant] == 'anc2c':
            anc2c.append(l_pc)
        elif l_n[variant] == 'anc2d':
            anc2d.append(l_pc)
        elif l_n[variant] == 'anc2e':
            anc2e.append(l_pc)
        elif l_n[variant] == 'anc2f':
            anc2f.append(l_pc)
        elif l_n[variant] == 'anc2g':
            anc2g.append(l_pc)
        elif l_n[variant] == 'anc2h':
            anc2h.append(l_pc)
        elif l_n[variant] == 'anc2i':
            anc2i.append(l_pc)
        elif l_n[variant] == 'anc2j':
            anc2j.append(l_pc)
        elif l_n[variant] == 'anc2k':
            anc2k.append(l_pc)
        elif l_n[variant] == 'anc2l':
            anc2l.append(l_pc)
        elif l_n[variant] == 'anc2m':
            anc2m.append(l_pc)
        elif l_n[variant] == 'anc2n':
            anc2n.append(l_pc)
        elif l_n[variant] == 'anc3':
            anc3.append(l_pc)
        elif l_n[variant] == 'anc4':
            anc4.append(l_pc)
        elif l_n[variant] == 'xcl1':
            xcl1.append(l_pc)



# In order to view the pLDDT scores of just the dimer files
#print(dimer_file_names)



colors=["red","gold","black"]
labels = ['', '', '']


# Creating Multiple Subplots for Pie Charts
fig, axes = plt.subplots(nrows=len(xcl1), ncols=19, figsize=(20, 3.5))

for i,j in zip(range(0,len(anc0)),anc0):
    axes[i,0].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2)),anc2):
    axes[i,1].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2a)),anc2a):
    axes[i,2].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2b)),anc2b):
    axes[i,3].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2c)),anc2c):
    axes[i,4].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2d)),anc2d):
    axes[i,5].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2d)),anc2e):
    axes[i,6].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2f)),anc2f):
    axes[i,7].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2g)),anc2g):
    axes[i,8].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2h)),anc2h):
    axes[i,9].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2i)),anc2i):
    axes[i,10].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2j)),anc2j):
    axes[i,11].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2k)),anc2k):
    axes[i,12].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2l)),anc2l):
    axes[i,13].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2m)),anc2m):
    axes[i,14].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc2n)),anc2n):
    axes[i,15].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc3)),anc3):
    axes[i,16].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(anc4)),anc4):
    axes[i,17].pie(x=j, labels=labels, colors=colors)

for i,j in zip(range(0,len(xcl1)),xcl1):
    axes[i,18].pie(x=j, labels=labels, colors=colors)



plt.show()

                

