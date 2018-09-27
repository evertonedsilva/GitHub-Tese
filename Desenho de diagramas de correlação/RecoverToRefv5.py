# -*- coding: utf-8 -*-
"""
Created on Thu May 10 13:25:42 2018

@author: Everton
"""
# libraries
#import pandas as pd
#import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import csv
import re


def returnSeq(fileSeq):
    seqFinal=''
    for eachLine in fileSeq:
        if not ">" in eachLine:        
            seqFinal+=eachLine.strip('\n')
    return seqFinal

#--------------------------------------------------
def EquivInRef(dic,ref,align,consensus):
    
    #print consensus 
    cons =returnSeq(consensus)
    
    
    ref =returnSeq(ref)
    align =returnSeq(align)
    
    alignungap=align.replace('-','')    
    alignungap = alignungap.upper()
    
    countRef = ref.find(alignungap)  
    
    cons=cons.upper()
    
    countAlign=cons.find(align.upper())
    
    begin = countAlign
    
        
    for aa in align: 
        countAlign+=1
        if aa != '-':
            dic.update({countAlign:{'aa':aa,'ref':countRef+1}})
            countRef+=1

    return (dic, [begin,begin+len(align)])

#------------------------------------------
dictRef={}

inputAlign = open('ArquivosRef/ConsensoSTMeTBDePleD_galperin.mfa','r')
ref1=open('ArquivosRef/ref1.txt','r')
align1=open('ArquivosRef/align1.txt','r')
dictRef= EquivInRef(dictRef,ref1,align1,inputAlign)
dom1= dictRef[1]
dictRef=dictRef[0]

#print len(dictRef)
inputAlign = open('ArquivosRef/ConsensoSTMeTBDePleD_galperin.mfa','r')
align2=open('ArquivosRef/align2.txt','r')
ref2=open('ArquivosRef/ref2.txt','r')
dictRef=EquivInRef(dictRef,ref2,align2,inputAlign)
dom2= dictRef[1]
dictRef=dictRef[0]
#print len(dictRef)

#print dom1
#print dom2

correlResult = {}

ID = []

color = {}

colordom1='green'
colordom2='darkorange'

with open('ArquivosRef/correlation_listv2.csv','r') as dataCorrel:
    fileCSV = csv.reader(dataCorrel)
    for row in fileCSV:        
        aa1 = re.findall(r'[A-Z]{1}', row[0])[0]
        num1=re.findall(r'[1-9].*', row[0])[0]
        num1=int(num1)
        if int(num1) in dictRef:
            aa1 = aa1+ str(dictRef[num1]['ref'])+'('+dictRef[num1]['aa']+')'
            #aa1 = aa1+ '('+str(num1)+')'
            
        else:
            aa1 = aa1+'('+str(num1)+')'
            
#adicionar a cor do nodo-----------------        
        if not aa1 in ID:
            ID.append(aa1)
            if dom1[0]<=num1<= dom1[1]:
                color.update({aa1:colordom1})
            elif dom2[0]<num1< dom2[1]:
                color.update({aa1:colordom2})
            else:
                color.update({aa1:'gray'})
       
            
        aa2 = re.findall(r'[A-Z]{1}', row[1])[0]
        num2=re.findall(r'[1-9].*', row[1])[0]
        num2=int(num2)
        if int(num2) in dictRef:
            aa2 = aa2+ str(dictRef[num2]['ref'])+'('+dictRef[num2]['aa']+')'
            #aa2 = aa2+'('+str(num2)+')'
            
            
        else:
            aa2 = aa2+'('+str(num2)+')'
        
             
        
        
#adicionar a cor do nodo ao dicionario----------------------------------------------
        if not aa2 in ID:
            ID.append(aa2)
            if dom1[0]<num2< dom1[1]:
                color.update({aa2:colordom1})
            elif dom2[0]<num2< dom2[1]:
                color.update({aa2:colordom2})
            else:
                color.update({aa2:'gray'})


#adicionar na lista de  pares correlacionados, considerando a volta (score passa a ser a media) 
        score = row[2]
        
        elem = (aa1,aa2)
        elemInv = (aa2,aa1)
       
        if not (elemInv in correlResult or elem in correlResult):
            correlResult.update({elem:{'corr':score}})
        
        else:
            score2 = correlResult[elemInv]['corr']
            
            correlResult[elemInv]['corr']=int((float(score)+float(score2))/2)
            
            


######.............................iniciando o desenho do rede -----------------------------

networkcorr = correlResult.keys()
#teste-----------------------------------------
networktest=[('L438(V)', 'F325(L)'), 
             ('A373(V)', 'F325(L)'), 
             ('K332(K)', 'D327(D)'), 
             ('D327(D)', 'Y439(Y)'), 
             ('D327(D)', 'H340(H)'), 
             ('Y439(Y)', 'L328(I)')] 
#networkcorr=networktest
#---------------------------------------------

G=nx.Graph()
G.add_edges_from(networkcorr)
# o pos eh opcional, mais interações=> mais proximos os nodos
# por experiencia: deixa sem pos
pos=nx.spring_layout(G,iterations=40)


# edges = conections---------------------------
scoresColor=[]
scores=[]
div=50
for corr in G.edges():
# por algum motivo essa merda fica trocando a ordem dos pares, por isso o 
# if e else seguinte:    
    if not corr in correlResult:
        corrEq = (corr[1],corr[0])
    else:
        corrEq = corr
    
    if int(correlResult[corrEq]['corr']) > 0:
        scoresColor.append('blue')

    else:
        scoresColor.append('red')
    
    scores.append(abs(float(correlResult[corrEq]['corr'])/div))
#-----------------------------------------------
# nodes = residues--------------------------------      
nodeColor=[]
for node in G.nodes():
    nodeColor.append(color[node])
#--------------------------------------------------


######......as proximas linha sao para testes (nao esquece de mudar em cima tambem)------
#networkcorr={('L438(V)', 'F325(L)'):20, ('A373(V)', 'F325(L)'):20, 
#             ('K332(K)', 'D327(D)'):1, ('D327(D)', 'Y439(Y)'):5, 
#            ('D327(D)', 'H340(H)'):5, ('Y439(Y)', 'L328(I)'):5}

#scores = []
#for corr in G.edges(): 
    #### por algum motivo essa merda fica trocando a ordem dos pares, por isso o 
###    if e else seguinte:
#     if not corr in networkcorr:
#        corrEq = (corr[1],corr[0])
#     else:
#        corrEq = corr       
        
#     scores.append(networkcorr[corrEq])


#nodes= ['H340(H)', 'A373(V)', 
#        'F325(L)', 'K332(K)', 
#        'D327(D)', 'Y439(Y)', 
#        'L328(I)', 'L438(V)'] 

#nodeColor=['green', 'green', 
#          'green', 'orange', 
#          'green', 'green', 
#          'purple', 'yellow']

#scoresColor=['gray', 'red', 
#             'green', 'blue', 
#             'gray', 'blue']
######......----------------------------------------------------------------------
    



#configuracoes e plotagem do grafico final----------------------------
    
# se tiver tudo muito junto, aumente o figsize ;D
    # melhor fazer tudo pequeno pro programa ir mais rapido
    
    
plt.figure(figsize=(30,30))
nx.draw_networkx(G, with_labels=True, 
        node_color=nodeColor, 
        node_size=100, 
        edge_color=scoresColor, 
        width=scores,
        font_size=2.2, 
        font_family='arial',
        font_weight='bold',
        linewidths=0.3)


ax = plt.gca() # to get the current axis
ax.collections[0].set_edgecolor("black")
plt.axis('off')
plt.savefig("FinalCorrGraph.svg", format="SVG")


# saída dos dados convertidos nas referências
Outputfile= open('correlacaoConv.txt','w')
Outputfile.write('aa1 \t aa2 \t log(p) \n')
outputOrder={'group1-1':[],'groupmix':[],'group2-2':[]}
for corr in G.edges():
# por algum motivo essa merda fica trocando a ordem dos pares, por isso o 
# if e else seguinte:    
    if not corr in correlResult:
        corrEq = (corr[1],corr[0])
    else:
        corrEq = corr
        
    aa1= corr[0]
    aa2 = corr[1]
    
    
    if color[aa1]!=color[aa2]:
        outputOrder['groupmix'].append(corrEq[0]+'\t'+corrEq[1]+'\t'+correlResult[corrEq]['corr']+'\n')
        
    else:       
        
        if color[aa1]==colordom1:
            outputOrder['group1-1'].append(str(corrEq[0])+'\t'+str(corrEq[1])+'\t'+str(correlResult[corrEq]['corr'])+'\n')
        else:
            outputOrder['group2-2'].append(corrEq[0]+'\t'+corrEq[1]+'\t'+correlResult[corrEq]['corr']+'\n')


#salve no arquivo de saida por grupos
for group in outputOrder:
    #apenas salve um marcador do grupo se ele nao for vazio
    if len(outputOrder[group])>0:
        Outputfile.write(group+'\n')
    #salve os elementos do grupo
    for ans in outputOrder[group]:
        Outputfile.write(ans)
  
    
Outputfile.close()
inputAlign.close()
ref1.close()
align1.close()
ref2.close()
align2.close()

#corrigir  variavel ID!!