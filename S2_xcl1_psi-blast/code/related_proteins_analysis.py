'''
Agglomerative clustering scatter plot for the PSI-blast hits for 1j9o
For S2
'''

# 60 pdbs total in folder 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.cluster import AgglomerativeClustering 



def create_matrix():
    m = np.zeros((60,60))

    names = []
    names_path = './S2_xcl1_psi-blast/psi-blast_hits/protein_names.txt'
    names_file = open(names_path, 'r')
    names_content = names_file.readlines()
    for names_line in names_content:
        names.append(names_line.strip())


    #### Note ###
    # Flip the i and the j to use the second TMscore
    j = 0
    for n in names:
        i = 0
        path = './S2_xcl1_psi-blast/psi-blast_hits/'+n+'_tmscores.txt'
        file = open(path,'r')

        content = file.readlines()

        for line in content:
            if i==j:
                    m[i][j] = np.nan
                    i+=1
            elif "(if normalized by length of Chain_1" in line:
                score = line.partition("= ")[2]
                score = score.partition(' (if normalized')[0]
                m[i][j] = score
                i+=1
        file.close() 
        j+=1   

    df = pd.DataFrame(m)
    df.columns = names
    df.to_csv('./S2_xcl1_psi-blast/psi-blast_hits/tmscore_matrix.csv')
    return df, names



def run_agglomerative(df, names):
    df = df.fillna(1)
    scale = StandardScaler()
    scale.fit(df)
    scaled_df = scale.transform(df)


    ac = AgglomerativeClustering(n_clusters = 2)
    clustered = ac.fit_predict(df)
    plt.scatter(df['1j9o'], df['2jp1'], edgecolors='black', c=clustered, cmap='autumn_r')
    plt.axhline(y=0.5, linestyle='--', c='gray')
    plt.axvline(x=0.5, linestyle='--', c='gray')

    # Using the first TM score
    plt.plot(0.31078,0.60120, label='2n54', c='black', marker='*', alpha=0.5, markersize=3) 
    plt.plot(0.68320,0.38434, label='1j8i', c='black', marker='*', alpha=0.5, markersize=3) 
    plt.plot(1,0.40167, label='1j9o', c='black', marker='*', alpha=0.5, markersize=3) 
    plt.plot(0.29073,1, label='2jp1', c='black', marker='*', alpha=0.5, markersize=3) 
    plt.plot(0.54925,0.37833, label='2hdm', c='black', marker='*', alpha=0.5, markersize=3) 
    x = [0.31078,0.68320,1,0.29073,0.54925]
    y=[0.60120,0.38434,0.40167,1,0.37833]



    label = ['2n54', '1j8i', '1j9o', '2jp1', '2hdm']
    for x,y, l in zip(x,y, label):
        plt.annotate(l, # this is the text
                    (x,y), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,10), # distance from text to points (x,y)
                    ha='center',
                    fontsize=12) # horizontal alignment can be left, right or center


    plt.tick_params(axis='x', labelsize=12)
    plt.tick_params(axis='y', labelsize=12)
    plt.xlabel('TM-score compared to 1j9o', fontsize=12)
    plt.ylabel('TM-score compared to 2jp1', fontsize=12)
    # To label every single point
    # for i, txt in enumerate(names):
    #     plt.annotate(txt, (df['1j9o'][i], df['2jp1'][i]))

    plt.show()

  



def main():
    df, names = create_matrix()
    run_agglomerative(df, names)


if __name__=="__main__":
    main()
    