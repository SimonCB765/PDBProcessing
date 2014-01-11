import gzip
import re

def main(mmCIFFile):
    """Parse an mmCIF file.

    Only parses a single file. Parses files that record only one PDB entry in them.

    @param mmCIFFile: The location of the gzipped mmCIF file to parse.
    @type mmCIFFile: string

    """

    readFile = gzip.open(mmCIFFile, 'r')
    mmCIFContent = readFile.read().decode('utf-8')
    readFile.close()

    fileChunks = re.split('(?<=\n)# \n', mmCIFContent)[:-1]

    entryID = fileChunks[1].split()[1]
    dataDictionary = {}
    entityIDs = []

    for i in fileChunks[1:]:
        # As the lines all end with newline characters, the final list element is always ''. [:-1] removes this.
        # It is necessary to do it this way, rather than splitting (i.split()), as this method ensures that all lines that should end with a ' ' do so.
        blockChunks = i.split('\n')[:-1]
        if blockChunks[0] == 'loop_':
            # A loop has been found, and therefore the block contains more than one data record.
            blockDictionary = {}
            # Record the name of the block.
            blockName = re.match('_[a-zA-Z0-9_]+', blockChunks[1])
            blockName = blockName.group(0)
            subBlockNameOrder = []  # Stored the subBlockNames in the order that they are encountered.
            currentNameOrderIndex = 0
            subBlockData = ''
            lookingAtLongDataReocrd = False
            longDataRecordDelimiter = ''
            secondLongDataReocrdDelimiter = False
            lookingAtMultilineDataReocrd = False
            dictSetUp = False
            for j in blockChunks[1:]:
                subBlockNameSearch = re.search('(?<=' + blockName + '.)[a-zA-Z0-9_]+', j)
                if subBlockNameSearch:
                    # If this is True, then the line contains something like _entity.id or _cell.angle_gamma_esd.
                    subBlockNameOrder.append(subBlockNameSearch.group(0))
                else:
                    # If this is True, then the line contains a data record.
                    if not dictSetUp:
                        blockDictionary = dict([(k, []) for k in subBlockNameOrder])
                        dictSetUp = True
                    previousCharacter = ''
                    for k in j:
                        # Go through every character on the line.
                        if lookingAtMultilineDataReocrd:
                            if previousCharacter == '' and k == ';':
                                # If this is True, then the beginning of a new line has been found, and the first character on the new line is a ;.
                                # This means that the end of a multiline data record has been found.
                                blockDictionary[subBlockNameOrder[currentNameOrderIndex]].append(subBlockData)
                                lookingAtMultilineDataReocrd = False
                                currentNameOrderIndex = (currentNameOrderIndex + 1) % len(subBlockNameOrder)
                                subBlockData = ''
                            else:
                                # If this is True, then the current character is in the middle of a multiline data record.
                                subBlockData += k
                            previousCharacter = k
                        elif lookingAtLongDataReocrd:
                            if k == longDataRecordDelimiter:
                                # If this is True, then the end of a '' or "" delimited data record has been found.
                                subBlockData += k
                                secondLongDataReocrdDelimiter = True
                            elif previousCharacter == longDataRecordDelimiter and k == ' ' and secondLongDataReocrdDelimiter:
                                lookingAtLongDataReocrd = False
                                longDataRecordDelimiter = ''
                                blockDictionary[subBlockNameOrder[currentNameOrderIndex]].append(subBlockData[:-1])
                                currentNameOrderIndex = (currentNameOrderIndex + 1) % len(subBlockNameOrder)
                                subBlockData = ''
                            else:
                                # If this is True, then the current character is in the middle of a long data record.
                                subBlockData += k
                            previousCharacter = k
                        else:
                            if previousCharacter == '' and k == ';':
                                # If this is True, then the beginning of a multiline data record has been found.
                                lookingAtMultilineDataReocrd = True
                                previousCharacter = k
                            elif k == '\'' or k == '"':
                                # If this is True, then the beginning of a long data record has been found.
                                lookingAtLongDataReocrd = True
                                secondLongDataReocrdDelimiter = False
                                longDataRecordDelimiter = k
                                previousCharacter = k
                            elif k == ' ':
                                if subBlockData != '':
                                    # If this is True, then the end of the current data record has been found.
                                    blockDictionary[subBlockNameOrder[currentNameOrderIndex]].append(subBlockData)
                                    currentNameOrderIndex = (currentNameOrderIndex + 1) % len(subBlockNameOrder)
                                    subBlockData = ''
                                previousCharacter = k
                            else:
                                # If this is True, then the current character is in the middle of a data record.
                                subBlockData += k
                                previousCharacter = k
        else:
            # There is no loop, and therefore the block only contains one data record.
            blockDictionary = {}
            # Record the name of the block.
            blockName = re.match('_[a-zA-Z0-9_]+', blockChunks[0])
            blockName = blockName.group(0)
            lookingAtMultilineDataReocrd = False
            for j in blockChunks:
                subBlockNameSearch = re.search('(?<=' + blockName + '.)[a-zA-Z0-9_]+', j)
                if subBlockNameSearch:
                    # If this is True, then the line contains something like _entity.id or _cell.angle_gamma_esd.
                    subBlockName = subBlockNameSearch.group(0)
                    subBlockChunks = j.split(None, 1)
                    if len(subBlockChunks) > 1:
                        # If this is True, then there is data on the line along with the subBlockName.
                        subBlockData = subBlockChunks[1].strip()
                        if subBlockData[-1] == '\'' or subBlockData[-1] == '"':
                            # Strip off leading and trailing '' or "".
                            blockDictionary[subBlockName] = [subBlockData[1:-1]]
                        else:
                            blockDictionary[subBlockName] = [subBlockData]
                    else:
                        # There is no data on the line with the subBlockName. This means that the data is spread over
                        # multiple lines, and is flanked by semi-colons.
                        pass
                else:
                    # There is no subBlockName on the line. This means that the line contains some data to go along with
                    # the most recently found subBlockName.
                    if lookingAtMultilineDataReocrd:
                        if j.strip() == ';':
                            # If this is True, then the line ends a data record.
                            lookingAtMultilineDataReocrd = False
                            blockDictionary[subBlockName] = [subBlockData]
                        else:
                            subBlockData += j
                    elif (j[0] == '\'' and j[-2:] == '\' ') or (j[0] == '"' and j[-2:] == '" '):
                        # If this is true then the line is taken up with a long piece of data that is surrounded by '' or "".
                        blockDictionary[subBlockName] = [j.strip()[1:-1]]
                    elif j[0] == ';':
                        # If this is True, then the line starts a data record.
                        lookingAtMultilineDataReocrd = True
                        subBlockData = j[1:]
                    else:
                        # If this is True, then the line just contains a data record that is too long for one line, but does not need/use
                        # the '' or "" long data record delimiters.
                        blockDictionary[subBlockName] = [j.strip()]

        dataDictionary[blockName] = blockDictionary

    return dataDictionary