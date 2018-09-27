
import re

inputSeq = open('GGDEFalignFilt_LinkerAlign_EALalignFilt.fasta','r')
#--------------
def DictConvert(inputSeq):
  forbiden = ['','\n','#','\t']
  dictConst = {}


  for eachLine in inputSeq:

        if not eachLine in forbiden:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.strip(' ')
            

            if ">" in eachLine:
                #code = eachLine.strip('>')
                code = eachLine.split('_')[0]
                dictConst.update({code:{'seqName':eachLine,'seq':''}})

            else:
                eachLine.upper()
                dictConst[code]['seq'] += eachLine

  return dictConst

dictSeq = DictConvert(inputSeq)

#----------

searchIn = {'GGDEF':{'RxGGDEF':{209:[],210:[],211:['A','G','S'],212:['G'],213:['D','E'],214:['E'],215:['F','L']},
                       'D-D':{89:['D'],116:['D']}
                   },
            'sitio-I':{'RxxD':{147:['R'],169:['D']}
                   },
            
            'EAL':{'Mg1': {885:['E'],1011:['N'],1083:['E'],1139:['D']},
                   'Mg2': {1140:['D']},
                   'Catx':{1283:['E']},
                   'Est': {1086:['E']}
                   }
                 				 	 		 	                                       
            }
# lembrar de adicionar um virgula se adicionar mais um grupo

#se desejar resíduos a uma distancia fixa, sem contar gaps: #:(i,+num,j)
#caso de anticorrelacão um ou outro, os dois simultâneos ocorre rejeicao:'Nha':{7:(['C'],'and not',6,['V'])}



	 	 		 	                                       
            
               				 	 		 	                                       
                               
motifs=searchIn.keys()
#motifs= sorted(motifs)


#para cada sequencia e verificada todos os elementos de cada motivo dado em searchIn
for seqName in dictSeq:

    dictSeq[seqName].update({'domains':{}})


    for motif in motifs:

        dictSeq[seqName]['domains'].update({motif:{}})
        conservedMotif='yep'

        for element in searchIn[motif]:
            dictSeq[seqName]['domains'][motif].update({element:{'found':'','conserved':'yep'}})

            found = re.sub(r'[A-VY-Za-z-1-9+_]','*',element)
        

            positions = searchIn[motif][element].keys()


            for position in positions:

                if type(searchIn[motif][element][position])is tuple and type(searchIn[motif][element][position][1])is int:

                    distance=''
                    query1=''
                    query2=''
                    for specific in searchIn[motif][element][position]:


                        #(aa na posicao dada,distancia sem gaps do proximo res?duo,aa esperado)
                        if type(specific)is int:
                            distance=specific


                        elif query1=='':
                            query1=specific

                        elif query1!='' and distance!='':
                            query2 = specific
                            aa = dictSeq[seqName]['seq'][position-1]
                            #print aa

                            if found.find('*') >= 0:
                                found=found[:found.find('*')] + aa + found[found.find('*')+1:]
                            else:
                                found=found+aa

                            #caso nao encontre um dos aminoacidos criticos
                            if not aa.upper() in query1:
                            #caso tenha sido fornecido que naum h? aminoacido critico
                                if not searchIn[motif][element][position]==[]:
                                    dictSeq[seqName]['domains'][motif][element]['conserved']='nope'
                                    conservedMotif = 'nope'

                            i=0
                            aa='-'

                            distance+=1

                            while aa=='-':
                                i+=1
                                if dictSeq[seqName]['seq'][position-1+i]!='-':
                                    #print dictSeq[seqName]['seq'][position-1+i]
                                    distance-=1
                                    if distance==0:
                                        aa = dictSeq[seqName]['seq'][position-1+i]
                                        #print aa


                            if found.find('*') >= 0:
                                found=found[:found.find('*')] + aa + found[found.find('*')+1:]
                            else:
                                found=found+aa


                            #caso nao encontre um dos aminoacidos criticos

                            if not aa.upper() in query2:
                            #caso tenha sido fornecido que naum ha aminoacido critico
                                if not searchIn[motif][element][position]==[]:
                                    conservedMotif  = 'nope'
                                    dictSeq[seqName]['domains'][motif][element]['conserved']='nope'

                            #print seqName
                            #print found

                            #reset dos par?mentros
                            query1=query2
                            query2=''
                            distance=''

                elif type(searchIn[motif][element][position])is tuple and type(searchIn[motif][element][position][1])is str:
                    
                    aa = dictSeq[seqName]['seq'][position-1]
                    aa2= dictSeq[seqName]['seq'][searchIn[motif][element][position][2]-1]
                    
                    
                    

                    #montagem da sequencia encontrada na regiao do motivo
                    if found.find('*') >= 0:
                        found=found[:found.find('*')] + aa + found[found.find('*')+1:]
                    else:
                        found+=aa
                    
                    #---------------
                    if found.find('*') >= 0:
                        found=found[:found.find('*')] + aa2 + found[found.find('*')+1:]
                    else:
                        found+=aa2
                    
                       
                        
                    if searchIn[motif][element][position][1]=='and not':
                        
                        # Se as duas posicoes tiverem presentes ou as duas tiverem nao estiverem a resposta é a rejeicao.
                        if (aa in searchIn[motif][element][position][0])==(aa2 in searchIn[motif][element][position][3]):
                            conservedMotif  = 'nope'
                            dictSeq[seqName]['domains'][motif][element]['conserved']='nope'
                            
                            
                                  
                    
                    
                    
                #casos normais em que a verificacao eh posicao por posicao (sem tuplas)
                else:
                    aa = dictSeq[seqName]['seq'][position-1]

                    #montagem da sequencia encontrada na regiao do motivo
                    if found.find('*') >= 0:
                        found=found[:found.find('*')] + aa + found[found.find('*')+1:]
                    else:
                        found=found + aa

                    if '*' in ''.join(searchIn[motif][element][position]):
                        queryObl=searchIn[motif][element][position]
                        queryOpt=[]

                        for Opt in searchIn[motif][element][position]:
                            if '*' in Opt:
                                queryObl.pop(queryObl.index(Opt))
                                Opt=Opt.strip('*')
                                queryOpt+=[Opt]

                        #caso nao encontre um dos aminoacidos criticos
                        if not aa.upper() in queryObl and not searchIn[motif][element][position]==[]:
                            if not aa.upper() in queryOpt:
                                conservedMotif  = 'nope'
                                dictSeq[seqName]['domains'][motif][element]['conserved']='nope'


                    else:
                        #caso nao encontre um dos aminoacidos criticos
                        if not aa.upper() in searchIn[motif][element][position] and not searchIn[motif][element][position]==[]:
                            conservedMotif  = 'nope'
                            dictSeq[seqName]['domains'][motif][element]['conserved']='nope'



            dictSeq[seqName]['domains'][motif][element]['found']=found

        #print  seqName
        dictSeq[seqName]['domains'][motif]['conserved']= conservedMotif
        #print  dictSeq[seqName]['domains']
        #print "\n"
        #put motif and if is conserved or not in Dict or something



#print motifs

motifKeys=[]
for motif in motifs:
    if motifKeys==[]:
        motifKeys=[motif+' conser|',motif+' nonconser|']

    else:
        new=[]
        for each in motifKeys:
            new.append(motif+' conser|'+each)
            new.append(motif+' nonconser|'+each)

        motifKeys=new





FinalDict={}

for item in motifKeys:
    elementKeys=[]
    for motif in motifs:
        
        for element in searchIn[motif].keys():
        
            if not motif + ' conser' in item:
                if elementKeys==[]:
                    elementKeys=[element+' conser|',element+' nonconser|']

                else:
                    new=[]
                    for each in elementKeys:
                        new.append(element+' conser|'+each)
                        new.append(element+' nonconser|'+each)
                    elementKeys=new
            else:
                if elementKeys==[]:
                    elementKeys=[element+' conser|']

                else:
                    new=[]
                    for each in elementKeys:
                        new.append(element+' conser|'+each)
                    elementKeys=new

        if motif + ' nonconser' in item:
            new =[]
            for verify in elementKeys:
                if not ' conser|'.join(searchIn[motif].keys()[::-1])+ ' conser' in verify:
                    new+=[verify]
            elementKeys = new

    FinalDict.update({item:{'num':0,'elements':{}}})
    for key in elementKeys:
        FinalDict[item]['elements'].update({key:{'found':{},'num':0}})



for seqName in dictSeq:
    motifkey=''
    elementkey=''
    found=[]
    for motif in dictSeq[seqName]['domains']:
        if dictSeq[seqName]['domains'][motif]['conserved']=='yep':
            motifkey=motif + ' conser|'+motifkey

        else:
            motifkey=motif + ' nonconser|'+motifkey

        for element in searchIn[motif].keys():
            found=[dictSeq[seqName]['domains'][motif][element]['found']]+found
            if dictSeq[seqName]['domains'][motif][element]['conserved']=='yep':
                elementkey=element + ' conser|'+elementkey

            else:
                elementkey=element + ' nonconser|'+elementkey
    #print found
    #print elementkey
    #print motifkey
    #print FinalDict
    FinalDict[motifkey]['num']+=1

    FinalDict[motifkey]['elements'][elementkey]['num']+=1

    found='...'.join(found)

    if not found in FinalDict[motifkey]['elements'][elementkey]['found']:
        FinalDict[motifkey]['elements'][elementkey]['found'].update({found:{'num':1,'codes':[]}})
    else:
        FinalDict[motifkey]['elements'][elementkey]['found'][found]['num']+=1

    code=seqName.split('/')[0]
    code=code.strip('>')
    
    #code=seqName.split('|')[2]
    #code=code.split(' ')[0]
    #print code


    FinalDict[motifkey]['elements'][elementkey]['found'][found]['codes'].append(code)
    #print FinalDict[motifkey]['elements'][elementkey]['found'][found]['codes']

summary = open('resultados/result_summary.txt', 'w')
for layers in FinalDict:
    name=layers.split('|')
    name=' '.join(name)
    OutputFile = open('resultados/'+name+'_com '+str(FinalDict[layers]['num'])+' seqs.txt', 'w')
    summary.write ('\t*****' + layers + ' ('+str(FinalDict[layers]['num'])+' )\n')
    print ('*****' + layers + ' ('+str(FinalDict[layers]['num'])+' )\n')
    
    for layer in FinalDict[layers]['elements']:
        if FinalDict[layers]['elements'][layer]['num']>0:
        
            summary.write ('-' + layer + ' ('+ str(FinalDict[layers]['elements'][layer]['num'])+') \n')
            for sublayer in FinalDict[layers]['elements'][layer]['found']:
                if FinalDict[layers]['elements'][layer]['found'][sublayer]['num']>0:
                    summary.write ( 'motivo: '+sublayer + ' ('+ str(FinalDict[layers]['elements'][layer]['found'][sublayer]['num'])+'):\n')
                    for code in FinalDict[layers]['elements'][layer]['found'][sublayer]['codes']:
                        OutputFile.write(">"+code+"\n")
                        OutputFile.write(dictSeq['>'+code]['seq']+'\n')
                        
                    summary.write ( '>'+' >'.join(FinalDict[layers]['elements'][layer]['found'][sublayer]['codes']))
                    summary.write ( '\n')
            summary.write ( '\n')
    summary.write ( '\n')
    
    OutputFile.close()



summary.close()
inputSeq.close()