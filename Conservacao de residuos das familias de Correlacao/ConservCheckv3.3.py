# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 13:01:50 2018

@author: Everton
"""
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



#--------------
def DictConvert(inputSeq):
  forbiden = ['','\n','#','\t']
  dictConst = {}

  for eachLine in inputSeq:

        if not eachLine in forbiden:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.strip(' ')

            if ">" in eachLine:
                code = eachLine
                dictConst.update({code:{'seq':''}})

            else:
                eachLine.upper()
                dictConst[code]['seq'] += eachLine

  return dictConst
#-------------------------------------------------------


#------------------------------------------------------------

def PrintTable(nodes,dictSeq,dictfull,dictRef,Outputfile):
    alignLen=len(dictSeq)
    alignLenFull=len(dictfull)
    positions = nodes.keys()
    positions = sorted(positions)
    possibilities='-PGAVILFWYCMSTDEKRHQN'
    possibilities=list(possibilities)
    zeroCount=[0]*len(possibilities)
    zeroCount=dict(zip(possibilities,zeroCount))
    


    for position in positions:
        nodes[position]={'aa':nodes[position],'found':zeroCount.copy(),'foundFull':zeroCount.copy()}
        for code in dictSeq:
            aa = dictSeq[code]['seq'][position-1]
            aa = aa.upper()
            if aa in possibilities:
                nodes[position]['found'][aa]+=1
            #print nodes[position]
            
    for position in positions:
        for code in dictfull:
            aa = dictfull[code]['seq'][position-1]
            aa = aa.upper()
            if aa in possibilities:
                nodes[position]['foundFull'][aa]+=1

    allpositions='\t'
    stringAux ='\t'
    for node in positions:
        allpositions+= nodes[node]['aa']+ str(dictRef[node]['ref'])+'('+dictRef[node]['aa'] +')\t\t'
        stringAux+= 'Grupo\tDelta\t'
    Outputfile.write(allpositions)
    Outputfile.write('\n')
    Outputfile.write(stringAux)
    Outputfile.write('\n')
    print allpositions
    print stringAux

    for possib in possibilities:
        allresult=possib+'\t'
        for node in positions:
            result =float(nodes[node]['found'][possib])
            result = result*100/alignLen
            
            result2 =float(nodes[node]['foundFull'][possib])
            result2 = result2*100/alignLenFull
            result2 = result - result2
            
            result = round(result,2)
            result2 = round(result2,2)
            
            allresult += str(result)+'\t'
            allresult += str(result2)+'\t'
        Outputfile.write(allresult)
        Outputfile.write('\n')
        
        print allresult
#------------------------------------------------
      
#------------------------------------------------------------

def PrintTableCompac(nodes,dictSeq,dictfull,dictRef,Outputfile):
    alignLen=len(dictSeq)
    alignLenFull=len(dictfull)
    positions = nodes.keys()
    positions = sorted(positions)
    possibilities='-PGAVILFWYCMSTDEKRHQN'
    possibilities=list(possibilities)
    zeroCount=[0]*len(possibilities)
    zeroCount=dict(zip(possibilities,zeroCount))
    


    for position in positions:
        nodes[position]={'aa':nodes[position],'found':zeroCount.copy(),'foundFull':zeroCount.copy()}
        for code in dictSeq:
            aa = dictSeq[code]['seq'][position-1]
            aa = aa.upper()
            if aa in possibilities:
                nodes[position]['found'][aa]+=1
            #print nodes[position]
            
    for position in positions:
        for code in dictfull:
            aa = dictfull[code]['seq'][position-1]
            aa = aa.upper()
            if aa in possibilities:
                nodes[position]['foundFull'][aa]+=1

    allpositions='\t'
    stringAux ='\t'
    for node in positions:
        allpositions+= nodes[node]['aa']+ str(dictRef[node]['ref'])+'('+dictRef[node]['aa'] +')\t\t\t'
        stringAux+= 'Grupo\tDelta\tMaxDelta\t'
    Outputfile.write(allpositions)
    Outputfile.write('\n')
    Outputfile.write(stringAux)
    Outputfile.write('\n')
    print allpositions
    print stringAux
    
    allresult='\t'
    for node in positions:
        possib = nodes[node]['aa']
        result =float(nodes[node]['found'][possib])
        result = result*100/alignLen
            
        result2 =float(nodes[node]['foundFull'][possib]) 
        result2 = result2*100/alignLenFull
        result2 = result - result2 #delta
        
        result = round(result,2)
        result2 = round(result2,2)
            
        allresult += str(result)+'\t'
        allresult += str(result2)+'\t'
        
        result3 = 0
        aaMax = '-'
        for possib in possibilities:        
            result =float(nodes[node]['found'][possib])
            result = result*100/alignLen
            
            result2 =float(nodes[node]['foundFull'][possib])
            result2 = result2*100/alignLenFull
            result2 = result - result2            

            result2 = round(result2,2)
            
            if abs(result2) > abs(result3):
                result3=result2
                aaMax = possib
    
        allresult += str(result3)+'('+aaMax+')'+'\t'     
            
        
    Outputfile.write(allresult)
    Outputfile.write('\n')
    print allresult      
        
#-------------------------------------------------------
def MaxDeltaConserv(alignBase,dictSeq,dictfull,dictRef,Outputfile3):
    nodes = {}
    i=0
    for position in alignBase:
        nodes.update({i:position})
        i+=1
        
    alignLen=len(dictSeq)
    alignLenFull=len(dictfull)
    positions = nodes.keys()
    positions = sorted(positions)
    possibilities='-PGAVILFWYCMSTDEKRHQN'
    possibilities=list(possibilities)
    zeroCount=[0]*len(possibilities)
    zeroCount=dict(zip(possibilities,zeroCount))   
    #nodes: posicoes do alinhamento

    for position in positions:
        nodes[position]={'aa':nodes[position],'found':zeroCount.copy(),'foundFull':zeroCount.copy()}
        for code in dictSeq:
            aa = dictSeq[code]['seq'][position-1]
            aa = aa.upper()
            if aa in possibilities:
                nodes[position]['found'][aa]+=1
            #print nodes[position]
            
    for position in positions:
        for code in dictfull:
            aa = dictfull[code]['seq'][position-1]
            aa = aa.upper()
            if aa in possibilities:
                nodes[position]['foundFull'][aa]+=1
    
    
    dictMaxResult = {}
    for node in positions:
        for possib in possibilities:        
            result =float(nodes[node]['found'][possib])
            result = result*100/alignLen
            
            result2 =float(nodes[node]['foundFull'][possib])
            result2 = result2*100/alignLenFull
            
            result3 = result - result2
        
            
            result = round(result,2)
            result2 = round(result2,2)
            result3 = round(result3,2)
            
            
            if result > 0 :
                
                ref='-'
                aaRef='-'
                
                #muitos das posicoes/nodes nao possui correspondente nas referencias
                if node in dictRef:
                    ref = dictRef[node]['ref']
                    aaRef = dictRef[node]['aa']
            
             
                resultName = possib + str(ref)+'('+aaRef+')'+'['+str(node)+']'+':\t'+str(result2)+'->'+str(result)
            
            
                if not abs(result3) in dictMaxResult:
                    dictMaxResult[abs(result3)]={'aas':[resultName]}
                else:
                    dictMaxResult[abs(result3)]['aas'].append(resultName)
                
    n=20
    maximus = sorted(dictMaxResult)
    
    maximus = maximus[::-1]

    for eachMax in maximus:
        for aa in  dictMaxResult[eachMax]['aas']:
            if n > 0:
                print aa
                n=n-1
            
            Outputfile3.write(aa)
            Outputfile3.write('\n')
        
        
    
    
#-------------------------------------------------------

        
inputSeq = open('ArquivosRef\EAL nonconser GGDEF nonconser _com 311 seqs.txt','r')
dictSeq = DictConvert(inputSeq)
inputAlignFull = open('ArquivosRef\GGDEFalignFilt_LinkerAlign_EALalignFilt.fasta','r')
dictSeqFull = DictConvert(inputAlignFull)
alignRefGGDEF_EAL = open('ArquivosRef\ConsensoSTMeTBDePleD_galperin.mfa','r')
Communs = open('ArquivosRef\communitiesCorrigido.txt','r')



ComColection = {}
ComOrder = []
for row in Communs : 
    row=row.strip('\n')
    if 'nodes' in row:
        group = row
        ComColection.update({group:{}})
        ComOrder.append(group)
    else:
        if not row in ['','\n']:
            num=re.findall(r'[1-9].*', row)[0]
            num=int(num)
            aa = re.findall(r'[A-Z]{1}', row)[0]            
            ComColection[group].update({num:aa})
            

a = 1


if a==0:
    Outputfile= open('resultados\dadosCorrelacao.txt','w')        
    for comun in  ComOrder: 
        Outputfile.write(comun+'\n')
        PrintTable(ComColection[comun],dictSeq,dictSeqFull,dictRef,Outputfile)
        print '\n'
        
    Outputfile.close()





if a==1:
    Outputfile2= open('resultados\dadosCorrelacaoCompact.txt','w')
    for comun in  ComOrder: 
        Outputfile2.write(comun+'\n')
        PrintTableCompac(ComColection[comun],dictSeq,dictSeqFull,dictRef,Outputfile2)
        print '\n'
    Outputfile2.close()

if a==2:
    Outputfile3= open('resultados\dadosVarConserv.txt','w')
    alignRef = returnSeq(alignRefGGDEF_EAL)
    
    MaxDeltaConserv(alignRef,dictSeq,dictSeqFull,dictRef,Outputfile3)

    Outputfile3.close()
           

    
inputSeq.close()

