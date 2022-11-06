import pandas as pd
import numpy as np
import math

dfGenotipo = pd.read_table('/IFTM_HACKATON/26k.txt')
dfMapa = pd.read_table('/IFTM_HACKATON/mapa_26k.txt')

colunas = dfGenotipo.columns.values.tolist()

indicesColunas = []


for j in range(len(colunas)):
    n = len(colunas[j])
    #print(n)
    aux = colunas[j]
    #print(aux)
    if(len(aux)>3):
        #print(aux[n-2])
        #print(aux[n-1])
        
        if(colunas[j].__contains__('Name')):
            print(aux)
            y = 'SNP Name'
            dfGenotipo = dfGenotipo.rename(columns={aux:y})

for j in range(len(colunas)):
    n = len(colunas[j])
    #print(n)
    aux = colunas[j]
    #print(aux)
    if(len(aux)>3):
        #print(aux[n-2])
        #print(aux[n-1])
        if(colunas[j].__contains__('Sample')):
            y = 'Sample ID'
            dfGenotipo = dfGenotipo.rename(columns={aux:y})

colunas = dfGenotipo.columns.values.tolist()




dfGenotipo2 = dfGenotipo.drop_duplicates(subset='SNP Name')

if(len(colunas)==4):
    dfGenotipo2 = dfGenotipo2.rename(columns={colunas[2]:'alelo 1'})
    dfGenotipo2 = dfGenotipo2.rename(columns={colunas[3]:'alelo 2'})
    
    flag = 0

    if(dfMapa.shape[0] == dfGenotipo2.shape[0]):
        #print(dfGenotipo2.shape[0])
        print('quantidade valida')
        
        i = 0;
       
        while(i<dfGenotipo2.shape[0] and flag == 0):
            print(i)
            if(dfMapa.loc[i,'Name'] != dfGenotipo2.loc[i,'SNP Name']):
                
                flag = 1
            i = i + 1
        

        if(flag==1):
            flag = 0
            print('Nomes invalidos')
            while(i<dfGenotipo2.shape[0] and flag == 0):
                print(i)
                if(dfMapa.loc[i,'Name'] != dfGenotipo2.loc[(dfGenotipo2.shape[0]-1)-i,'SNP Name']):
                    
                    flag = 1
                i = i + 1
            
        else:
            print('Nomes validados')
        
        if(flag == 1):
            print('Nomes invalidos')
        else:
            print('Nomes validados')
            
    
else:
    
    flag = 0

    if(dfMapa.shape[0] == dfGenotipo2.shape[0]):
        #print(dfGenotipo2.shape[0])
        print('quantidade valida')
        
        i = 0;
        
        while(i<dfGenotipo2.shape[0] and flag == 0):
            if(dfMapa.loc[i,'Name'] != dfGenotipo2.loc[i,'SNP Name']):
                flag = 1
            i = i + 1
        
        if(flag==1):
            print('Nomes invalidos')
            
        else:
            print('Nomes validados')


    print(colunas)

    for j in range(len(colunas)):
        n = len(colunas[j])
        #print(n)
        aux = colunas[j]
        #print(aux)
        if(len(aux)>3):
            #print(aux[n-2])
            #print(aux[n-1])
            if(aux[n-2]=='A' and aux[n-1]=='B'):
                indicesColunas.append(j)


    print(indicesColunas)
    




SNP_Name = np.asarray(dfGenotipo2['SNP Name'].str.replace('-','.'))
Sample_ID = np.asarray(dfGenotipo2['Sample ID'])
if(len(colunas)>4):
    alelo1 = np.asarray(dfGenotipo2[colunas[indicesColunas[0]]])
    alelo2 = np.asarray(dfGenotipo2[colunas[indicesColunas[1]]])
else:
    alelo1 = np.asarray(dfGenotipo2['alelo 1'])
    alelo4 = np.asarray(dfGenotipo2['alelo 2'])
gene = alelo1+alelo2

for j in range(len(gene)):
    cont= 0
    if(gene[j] == 'AA'):
        cont = cont + 1
        gene[j] = 0
    
    if(gene[j] == 'BB'):
        gene[j] = 2
        cont = cont + 1
    if(gene[j] == 'AB'):
        gene[j] = 1
        cont = cont + 1
    if(gene[j] == 'NRNR' or gene[j] == '--'):
        gene[j] = 5
        cont = cont + 1
    if(cont == 0):
        j= len(gene)
        print('Gene da linha ',j,' está incorreto.')
        
array = []
array.append(SNP_Name)
array.append(Sample_ID)
array.append(gene)

#array = np.transpose(array)

df3 = pd.DataFrame(data = array)
valor = math.floor(len(gene)/1000)

csv = df3.to_csv('genes'+str(valor)+'K.csv')

