import xml.etree.ElementTree as ET
import os


tree = ET.parse('KIR.xml')
root = tree.getroot()

print(root.tag)

#with open('FASTA_Format_Alleles.txt', 'w') as file:
os.system('rm FASTA_*')
for geneInfo in root:
    file = None
    geneFamily = None
    nucsequence = ""
    formattedSequence = ""
    borders = []
    exonsDone = []

    KIR2DL2_KIR2DL3_Flag = 0
    KIR2DS3_KIR2DS5_Flag = 0
    KIR3DL1_KIR3DS1_Flag = 0
    KIR2DL5_Flag = 0
    #file.write("Allele Name: " + geneInfo.get('name') + '\n')
    header = ">" + geneInfo.get('name') + "\n"
    typeLocus = None
    for attributes in geneInfo:
        if str(attributes.tag).find('sequence') > -1:
            for sequence in attributes:
                #if str(sequence.tag).find('alignmentreference') > -1:
                    #file.write("Allele Name: " + sequence.get('allelename') + "\n")
                
                #file.write(sequence.tag + "\n")
                if str(sequence.tag).find('nucsequence') > -1:
                    nucsequence = sequence.text
                    #file.write("General Sequence: " + nucsequence + "\n")

                if str(sequence.tag).find('feature') > -1:
                    for information in sequence:
                        if str(information.tag).find('translation') == -1:
                            coordinateTypes = information.tag[information.tag.find('}') + 1:]
                            exonType = sequence.get('featuretype')
                            exonName = sequence.get('name')
                            if exonName not in exonsDone:
                                exonsDone.append(exonName)
                                #file.write(exonType + ": " + exonName + "\n")
                            sequenceCoordinatesStart = int(information.get('start'))
                            sequenceCoordinatesEnd = int(information.get('end'))
                            #file.write(coordinateTypes + ":\t\t")
                            if coordinateTypes.find('cDNA') > -1:
                                #file.write("\t")
                                x = 0
                            else:
                                borders.append(sequenceCoordinatesStart - 1)
                                borders.append(sequenceCoordinatesEnd)
                            #file.write(str(sequenceCoordinatesStart) + "---" + str(sequenceCoordinatesEnd) + "\n")
                                
            tempSequence = ""

            for i in range(0, len(exonsDone)):
                if exonsDone[i].find('Exon') > -1 and exonsDone[i].find('Pseudo') == -1: #2i and 2i + 1
                    # KIR2DS Series incorrectly labeled in KIR.xml. Psuedoexon 3 is labeled as an Exon.
                    if typeLocus.find('KIR2DS') > -1 and exonsDone[i].find('3') > -1:
                        print('KIR.xml incorrectly labeled.')
                    else: 
                        rangeStart = borders[2*i]
                        rangeEnd = borders[2*i + 1]
                        tempSequence += nucsequence[rangeStart:rangeEnd]

            count = 0
            
            while count+60 <= len(tempSequence):
                formattedSequence += tempSequence[count:count+60] + "\n"
                count += 60
            formattedSequence += tempSequence[count:] + '\n'
            file.write(header + formattedSequence)
            
            if KIR2DL2_KIR2DL3_Flag == 1:
                with open('FASTA_KIR2DL2_KIR2DL3.fasta', 'a') as joint_file:
                    joint_file.write(header + formattedSequence)
                KIR2DL2_KIR2DL3_Flag = 0
            
            if KIR2DS3_KIR2DS5_Flag == 1:
                with open('FASTA_KIR2DS3_KIR2DS5.fasta', 'a') as joint_file:
                    joint_file.write(header + formattedSequence)
                KIR2DS3_KIR2DS5_Flag = 0

            if KIR3DL1_KIR3DS1_Flag == 1:
                with open('FASTA_KIR3DL1_KIR3DS1.fasta', 'a') as joint_file:
                    joint_file.write(header + formattedSequence)
                KIR3DL1_KIR3DS1_Flag = 0
            
            if KIR2DL5_Flag == 1:
                with open('FASTA_KIR2DL5.fasta', 'a') as joint_file:
                    joint_file.write(header + formattedSequence)
                KIR2DL5_Flag = 0

            with open('FASTA_General.fasta', 'a') as general_file:
                general_file.write(header + formattedSequence)
            #file.write("EXPERIMENTAL SEQUENCE:\n" + formattedSequence + "\n\n")

        if str(attributes.tag).find('locus') > -1:
            file = open('FASTA_' + attributes.get('locusname') + '.fasta', 'a')
            typeLocus = attributes.get('locusname')            
            if attributes.get('locusname') == 'KIR2DL2' or attributes.get('locusname') == 'KIR2DL3':
                KIR2DL2_KIR2DL3_Flag = 1
            if attributes.get('locusname') == 'KIR2DS3' or attributes.get('locusname') == 'KIR2DS5':
                KIR2DS3_KIR2DS5_Flag = 1
            if attributes.get('locusname') == 'KIR3DL1' or attributes.get('locusname') == 'KIR3DS1':
                KIR3DL1_KIR3DS1_Flag = 1
            if attributes.get('locusname') == 'KIR2DL5A' or attributes.get('locusname') == 'KIR2DL5B':
                KIR2DL5_Flag = 1