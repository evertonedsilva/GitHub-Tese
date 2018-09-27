import urllib
import re

webpage = urllib.urlopen('http://www.ncbi.nlm.nih.gov/Complete_Genomes/c-di-GMP.html').read()

webpage = re.findall(r'^<td><A href=".*?\r\n</tr>',webpage,re.DOTALL|re.MULTILINE|re.IGNORECASE)

webpage = webpage
#for found in webpage:
#    print found
#    print '\n'

seqTotal = 0
FoundSeq = 0
FailSeq = 0
OutputFile = open('ggdef_eal.txt', 'w')
OutputFail = open('ggdef_eal_broken.txt', 'w')
for found in webpage:
    if found != []:
        found = found.split('\r\n')

        if found[4]== '<td><center>-</td>':
            organism = found[0].split('">')[1]
            organism = organism.split('<')[0]
            GGDEF_EAL = 0

            #print organism + " (" + str(GGDEF_EAL) + ') :'

        else:
            orgLink=found[4].split('<center><A href=')[1]
            orgLink=orgLink.split(">")[0].strip('"')
            #print orgLink

            GGDEF_EAL=found[4].split(">")[3]
            GGDEF_EAL=int (GGDEF_EAL.split('<')[0])


            organism=found[0].split('">')[1]
            organism=organism.split('<')[0]
            #print organism + ' : '+ GGDEF_EAL
            seqTotal+= GGDEF_EAL



            orgIdGp=found[4].split('#')[1]
            #print found[4]
            orgIdGp=orgIdGp.split('">')[0]
            #print orgIdGp

            webpageOrg = urllib.urlopen('http://www.ncbi.nlm.nih.gov'+orgLink).read()
            #print 'http://www.ncbi.nlm.nih.gov'+orgLink
            #print orgIdGp
            #print webpageOrg


            webpageOrg = re.findall(r'<a name=\s??'+orgIdGp+'\s??>.*?<HR>',webpageOrg, re.DOTALL|re.IGNORECASE|re.MULTILINE)
            #print  orgIdGp, organism

            if webpageOrg == []:
                OutputFail.write('> Gaperin link code error ' + organism + ': \n link code:' + orgIdGp + '\n' )
                FailSeq += GGDEF_EAL
                #print '>Nao foi possivel acessar os links do organism ' + organism + ' usando o codigo ' + orgIdGp
                #print link

            else:

                webpageOrg = webpageOrg[0]
                Domain = ''
                GotAll = 1

                if not 'GGDEF-EAL' in webpageOrg:
                    Domain = 'GGDEF-EAL'
                    GotAll = 0

                print '\n' + organism + " (" + str(GGDEF_EAL) + ') :'
                print '- Link: http://www.ncbi.nlm.nih.gov'+orgLink
                print '- Links Uniprot com problema:'

                webpageOrg = webpageOrg.split('\n')
                #print webpageOrg

                GGDEF_EAL_page1 = GGDEF_EAL
                SeqCount = 0
                CountFail = 0

                for line in webpageOrg:
                    if '<a name=' in line or '<A name' in line:
                        pass
                    elif 'only' in line or 'HD-GYP' in line:
                        Domain = 'One Domain or HD-GYP'
                    elif 'GGDEF-EAL' in line and Domain!='GGDEF-EAL':
                        #isso faz com que a seja lida apenas uma vez as linhas da webpage a partir que encontra uma linha GGDEF-EAL
                        Domain = 'GGDEF-EAL'


                    elif not '<HR>' in line and Domain=='GGDEF-EAL':
                        GGDEF_EAL-=1
                        #print GGDEF_EAL
                        if GGDEF_EAL >= 0 or GotAll ==1:
                            #se GotAll == 1 quer dizer que ha o identificador GGDEF-EAL na linha entaum pegue todas as seqs mesmo se ha divergencia no numero de GGDEF-EAL
                            #print line
                            link = re.split(r'href="(?i)', line, re.DOTALL | re.IGNORECASE | re.MULTILINE)
                            #print link
                            #print organism
                            link = link[4].split('"')[0]
                            link = link.strip(' ')
                            link = link+'.fasta'

                            seq = urllib.urlopen(link).read()


                            if seq =='':
                                line = line.split(' ')[0]
                                OutputFail.write('> Uniprot link empity -' +organism + line + ': \n' + link + '\n')
                                FailSeq += 1
                                CountFail +=1
                                print link


                                #print '>A seq do organismo ' + organism + ' com codigo original ' + line + ' possui problema no link do Uniprot'
                                #print link

                            else:
                                OutputFile.write(seq)
                                FoundSeq += 1
                                SeqCount +=1
                                #print SeqCount

                        else:
                            print '- Este link foi descartado:'
                            print link

                if GGDEF_EAL != 0:
                    print '- Existe uma incoerencia entre numero de seqs indicadas com as recuperadas: '
                    print "Seqs da pag. principal = " + str (GGDEF_EAL_page1)
                    print 'Seqs recuperadas       = ' + str(SeqCount)
                    print 'Seqs com erro          = ' + str (CountFail)
                    print 'Links encontrados      = ' + str(GGDEF_EAL_page1 - GGDEF_EAL)



print 'Total de seqs da pagina principal = ' + str(seqTotal)
print 'Total de seqs as recuperadas      = ' + str(FoundSeq)
print 'Total de seqs com erro  = ' + str (FailSeq) + '\n'


OutputFile.close()
OutputFail.close()


