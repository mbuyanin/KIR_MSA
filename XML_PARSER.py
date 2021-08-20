import xml.etree.ElementTree as ET

tree = ET.parse('KIR.xml')
root = tree.getroot()

print(root.tag)
with open('PARSED_FILE.txt', 'w') as file:
    for geneInfo in root:
        nucsequence = ""
        borders = []
        exonsDone = []
        file.write("Allele Name: " + geneInfo.get('name') + '\n')
        for attributes in geneInfo:
            if str(attributes.tag).find('sequence') > -1:
                for sequence in attributes:
                    #if str(sequence.tag).find('alignmentreference') > -1:
                        #file.write("Allele Name: " + sequence.get('allelename') + "\n")
                    
                    #file.write(sequence.tag + "\n")
                    if str(sequence.tag).find('nucsequence') > -1:
                        nucsequence = sequence.text
                        file.write("General Sequence: " + nucsequence + "\n")

                    if str(sequence.tag).find('feature') > -1:
                        for information in sequence:
                            if str(information.tag).find('translation') == -1:
                                coordinateTypes = information.tag[information.tag.find('}') + 1:]
                                exonType = sequence.get('featuretype')
                                exonName = sequence.get('name')
                                if exonName not in exonsDone:
                                    exonsDone.append(exonName)
                                    file.write(exonType + ": " + exonName + "\n")
                                sequenceCoordinatesStart = int(information.get('start'))
                                sequenceCoordinatesEnd = int(information.get('end'))
                                file.write(coordinateTypes + ":\t\t")
                                if coordinateTypes.find('cDNA') > -1:
                                    file.write("\t")
                                else:
                                    borders.append(sequenceCoordinatesStart - 1)
                                    borders.append(sequenceCoordinatesEnd)
                                file.write(str(sequenceCoordinatesStart) + "---" + str(sequenceCoordinatesEnd) + "\n")
                                    
                tempSequence = ""
                for i in range(0, len(exonsDone)):
                    if exonsDone[i].find('Exon') > -1: #2i and 2i + 1
                        rangeStart = borders[2*i]
                        rangeEnd = borders[2*i + 1]
                        tempSequence += nucsequence[rangeStart:rangeEnd]

                count = 0
                formattedSequence = ""
                while count+60 <= len(tempSequence):
                    formattedSequence += tempSequence[count:count+60] + "\n"
                    count += 60
                formattedSequence += tempSequence[count:]
                file.write("EXPERIMENTAL SEQUENCE:\n" + formattedSequence + "\n\n")