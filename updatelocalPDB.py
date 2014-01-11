import os
import shutil
import sqlite3
import subprocess
import sys
import time

import parsePDBmmCIF
import processPSIoutput

def main(mmCIFDir, parsedPDB, blastExecutables, daysSinceLastUpdate=7):
    """
    """

    # Define output files.
    if not os.path.exists(parsedPDB):
        os.mkdir(parsedPDB)
    fileAllPDBEntries = parsedPDB + '/AllPDBEntries.txt'
    fileChainType = parsedPDB + '/ChainType.txt'
    fileProteinInformation = parsedPDB + '/ProteinInformation.txt'
    fileRepresentative = parsedPDB + '/Representative.txt'
    fileSimilarity = parsedPDB + '/Similarity.txt'

    # Determine the mmCIF files that need parsing.
    secondsInADay = 86400
    currentTime = time.time()
    mmCIFDirs = os.listdir(mmCIFDir)
    toParse = []
    for i in mmCIFDirs:
        mmCIFFolder = mmCIFDir + '/' + i
        mmCIFFiles = os.listdir(mmCIFFolder)
        for j in mmCIFFiles:
            currentFile = mmCIFFolder + '/' + j
            lastModified = os.path.getmtime(currentFile)
            if currentTime - (secondsInADay * daysSinceLastUpdate) < lastModified:
                # If the file has been created/modified within the last timeframe then it needs parsing.
                toParse.append(currentFile)

    # Connect to the SQLite database
    conn = sqlite3.connect('Leaf.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS ChainTypes (chain TEXT PRIMARY KEY, chainType TEXT)''')
    c.execute('''CREATE TABLE IF NOT EXISTS ProteinInformation (chain TEXT PRIMARY KEY, entry TEXT, experimentType TEXT, resolution REAL, rValueObs REAL,
                 rValueFree REAL, alphaCarbonOnly INTEGER, description TEXT, dbName TEXT, dbCode TEXT, organism TEXT, sequence TEXT)''')
    c.execute('''CREATE TABLE IF NOT EXISTS Similarities (chainA TEXT, entryA TEXT, chainB TEXT, entryB TEXT, similarity REAL, matchLength INTEGER,
                 PRIMARY KEY (chainA, chainB))''')
    conn.commit()

    typeDict = {}
    proteinDict = {}
    newlyParsedChains = set([])
    for i in toParse:
        # Parse mmCIF file.
        entryID, entityRecords, experimentalType, resolution, rFactorObs, rFactorFree = parsePDBmmCIF.main(i)

        # Chain type information.
        for j in entityRecords.keys():
            if 'type' in entityRecords[j]:
                chains = [entry + chain for entry in entryID for chain in entityRecords[j]['chains']]
                newlyParsedChains |= set(chains)
                type = entityRecords[j]['type'].strip()

                if type == 'Protein':
                    dbCode = entityRecords[j]['dbCode'].strip() if 'dbCode' in entityRecords[j] else ''
                    if dbCode in ['?', '.']:
                        dbCode = ''
                    dbName = entityRecords[j]['dbName'].strip() if 'dbName' in entityRecords[j] else ''
                    if dbName in ['?', '.']:
                        dbName = ''
                    description = entityRecords[j]['description'].strip() if 'description' in entityRecords[j] else ''
                    if description in ['?', '.']:
                        description = ''
                    onlyAlphaCarbon = entityRecords[j]['onlyAlphaCarbon']
                    scientificName = entityRecords[j]['scientificName'].strip() if 'scientificName' in entityRecords[j] else ''
                    if scientificName in ['?', '.']:
                        scientificName = ''
                    sequence = entityRecords[j]['sequence'].upper()
                    # Determine if over 50% of the amino acids are X, in which case PISCES marks the protein as being a nonprotein.
                    if sequence.count('X') / float(len(sequence)) >= 0.5:
                        type = 'NonProtein'
                    else:
                        for k in chains:
                            proteinDict[k] = {'dbCode' : dbCode, 'dbName' : dbName, 'description' : description, 'entry' : entryID[0],
                                             'experimentalType' : experimentalType, 'onlyAlphaCarbon' : onlyAlphaCarbon, 'resolution' : resolution,
                                             'rFactorFree' : rFactorFree, 'rFactorObs' : rFactorObs, 'scientificName' : scientificName,
                                             'sequence' : sequence}

                # Record the type of the entity.
                for k in chains:
                    typeDict[k] = type

    # Determine different subsets of the chains.
    c.execute('SELECT chain FROM ChainTypes')
    existingChains = set([i[0] for i in c.fetchall()])
    addedChains = newlyParsedChains - existingChains
    updatedChains = newlyParsedChains & existingChains
    deletedChains = existingChains - newlyParsedChains
    c.execute('SELECT chain FROM ChainTypes WHERE chainType=\'Protein\'')
    existingProteinChains = set([i[0] for i in c.fetchall()])

    ##########################
    # Update altered chains. #
    ##########################
    # Only care about updates when the sequence has changed (as this is the only thing that affects the similarity).
    chainsWithSimilaritiesToAlter = set([])
    for i in updatedChains:
        c.execute('SELECT sequence FROM ProteinInformation WHERE chain=?', (i,))
        result = c.fetchone()
        if result and result[0] == proteinDict[i]['sequence']:
            # If the sequences are the same, then just update the info in the record.
            c.execute('DELETE FROM ProteinInformation WHERE chain=?', (i,))
            tupleToAdd = (i, proteinDict[i]['entry'], proteinDict[i]['experimentalType'], proteinDict[i]['resolution'], proteinDict[i]['rFactorObs'],
                          proteinDict[i]['rFactorFree'], 1 if proteinDict[i]['onlyAlphaCarbon'] else 0, proteinDict[i]['description'], proteinDict[i]['dbName'],
                          proteinDict[i]['dbCode'], proteinDict[i]['scientificName'], proteinDict[i]['sequence'])
            c.execute('INSERT INTO ProteinInformation VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', tupleToAdd)
            conn.commit()
        elif not i in existingProteinChains:
            # If the chain is not an existing protein chain, then there is nothing to update.
            pass
        else:
            # Treat the chain as if it was just added (and also delete its information first).
            addedChains |= set([i])
            deletedChains |= set([i])

    ######################
    # Delete old chains. #
    ######################
    c.executemany('DELETE FROM ChainTypes WHERE chain=?', [(i,) for i in deletedChains])
    c.executemany('DELETE FROM ProteinInformation WHERE chain=?', [(i,) for i in deletedChains])
    c.executemany('DELETE FROM Similarities WHERE chainA=?', [(i,) for i in deletedChains])
    c.executemany('DELETE FROM Similarities WHERE chainB=?', [(i,) for i in deletedChains])
    conn.commit()

    ###################
    # Add new chains. #
    ###################
    # Add the type information for the new chains.
    typeTuplesToAdd = [(i, typeDict[i]) for i in addedChains]
    c.executemany('INSERT INTO ChainTypes VALUES (?, ?)', typeTuplesToAdd)

    # Check if any of the new protein chains are identical to ones already recorded.
    addedProteinChains = set([i for i in addedChains if typeDict[i] == 'Protein'])
    identicalAddedProteinChains = set([])
    for i in addedProteinChains:
        c.execute('SELECT chain FROM ProteinInformation WHERE sequence=?', (proteinDict[i]['sequence'],))
        result = c.fetchone()
        if result:
            # If the chain has an identical sequence to a chain that is already in the database, then copy its similarities.
            identicalAddedProteinChains.add(i)
            identicalChain = result[0]
            c.execute('SELECT * FROM Similarities WHERE chainA=?', (identicalChain,))
            chainASimilarities = [j[0] for j in c.fetchall()]
            for j in chainASimilarities:
                j[0] = i
            c.execute('SELECT * FROM Similarities WHERE chainB=?', (identicalChain,))
            chainBSimilarities = [j[0] for j in c.fetchall()]
            for j in chainBSimilarities:
                j[2] = i
            allChainSimilarities = set([])
            allChainSimilarities |= set(chainASimilarities)
            allChainSimilarities |= set(chainBSimilarities)
            c.executemany('INSERT INTO Similarities VALUES (?, ?, ?, ?, ?, ?)', allChainSimilarities)

    # Add the protein information for the new chains.
    protInfoTuplesToAdd = [(i, proteinDict[i]['entry'], proteinDict[i]['experimentalType'], proteinDict[i]['resolution'], proteinDict[i]['rFactorObs'],
                           proteinDict[i]['rFactorFree'], 1 if proteinDict[i]['onlyAlphaCarbon'] else 0, proteinDict[i]['description'], proteinDict[i]['dbName'],
                           proteinDict[i]['dbCode'], proteinDict[i]['scientificName'], proteinDict[i]['sequence']) for i in addedProteinChains]
    c.executemany('INSERT INTO ProteinInformation VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', protInfoTuplesToAdd)
    conn.commit()

    # As the identical ones already have their similarities recorded, they count as existing chains rather than newly added ones.
    addedProteinChains -= identicalAddedProteinChains
    existingProteinChains |= identicalAddedProteinChains

    #############################################
    # Determine sequences for and run BLASTing. #
    #############################################
    existingChainsFASTA = parsedPDB + '/ExistingChains.fasta'
    addedChainsFASTA = parsedPDB + '/AddedChains.fasta'
    if addedProteinChains:
        # Determine the sequences for the existing protein chains.
        c.execute('SELECT chain, sequence FROM ProteinInformation WHERE chain IN (\'' + '\',\''.join(existingProteinChains) + '\')')
        existingProteinsResult = c.fetchall()

        # Take only the unique existing chains.
        representativesExisting = {}
        for i in existingProteinsResult:
            chain = i[0]
            sequence = i[1]
            if sequence in representativesExisting:
                representativesExisting[sequence].append(chain)
            else:
                representativesExisting[sequence] = [chain]

        # Generate FASTA file of the existing protein chains.
        writeExisting = open(existingChainsFASTA, 'w')
        for seq in representativesExisting:
            chains = ','.join(representativesExisting[seq])
            writeExisting.write('>' + chains + '\n' + seq + '\n')
        writeExisting.close()

        # Determine the sequences for the added non-identical protein chains.
        c.execute('SELECT chain, sequence FROM ProteinInformation WHERE chain IN (\'' + '\',\''.join(addedProteinChains) + '\')')
        addedProteinsResult = c.fetchall()

        # Take only the unique added chains.
        representativesAdded = {}
        for i in addedProteinsResult:
            chain = i[0]
            sequence = i[1]
            if sequence in representativesAdded:
                representativesAdded[sequence].append(chain)
            else:
                representativesAdded[sequence] = [chain]

        # Generate FASTA file of added protein chains.
        writeAdded = open(addedChainsFASTA, 'w')
        for seq in representativesAdded:
            chains = ','.join(representativesAdded[seq])
            writeAdded.write('>' + chains + '\n' + seq + '\n')
        writeAdded.close()

        # Generate BLAST databases.
        existingDatabaseDir = parsedPDB + '/ExistingDatabase'
        os.mkdir(existingDatabaseDir)
        os.mkdir(existingDatabaseDir + '/TempDB')
        makeDBArgs = [blastExecutables + '/makeblastdb', '-in', existingChainsFASTA, '-out', existingDatabaseDir + '/TempDB', '-dbtype', 'prot']
        subprocess.call(makeDBArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        addedDatabaseDir = parsedPDB + '/AddedDatabase'
        os.mkdir(addedDatabaseDir)
        os.mkdir(addedDatabaseDir + '/TempDB')
        makeDBArgs = [blastExecutables + '/makeblastdb', '-in', addedChainsFASTA, '-out', addedDatabaseDir + '/TempDB', '-dbtype', 'prot']
        subprocess.call(makeDBArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # BLAST added against existing and vice versa.
        if existingProteinChains:
            # If there are already some proteins in the database (will only fail when the database is first being created).
            # BLAST added against existing.
            addedAgainstExisting = parsedPDB + '/AddedAgainstExisting.txt'
            sequence_BLAST(addedChainsFASTA, addedAgainstExisting, existingDatabaseDir + '/TempDB', blastExecutables + '/psiblast.exe', 2)
            similaritiesAddedAgainstExisting = processPSIoutput.main(addedAgainstExisting)

            #BLAST existing against added.
            existingAgainstAdded = parsedPDB + '/ExistingAgainstAdded.txt'
            sequence_BLAST(existingChainsFASTA, existingAgainstAdded, addedDatabaseDir + '/TempDB', blastExecutables + '/psiblast.exe', 2)
            similaritiesExistingAgainstAdded = processPSIoutput.main(existingAgainstAdded)

            # Determine similarities to add.
            for i in similaritiesAddedAgainstExisting:
                query = i[0].split(',')
                hits = i[1].split(',')
                similarity = similaritiesAddedAgainstExisting[i]['Identity']
                matchLength = similaritiesAddedAgainstExisting[i]['MatchLength']
                pairs = [[i,j] for i in query for j in hits]
                if i in similaritiesExistingAgainstAdded:
                    # The query/hit pair is in both dicts, so determine which dict holds the greater similarity.
                    alternateSimilarity = similaritiesExistingAgainstAdded[i]['Identity']
                    if alternateSimilarity > similarity:
                        similarity = alternateSimilarity
                        matchLength = similaritiesExistingAgainstAdded[i]['MatchLength']
                similarityTuples = [(i[0], i[0][:-1], i[1], i[1][:-1], similarity, matchLength) for i in pairs]
                c.executemany('INSERT INTO Similarities VALUES (?, ?, ?, ?, ?, ?)', similarityTuples)
                conn.commit()
            os.remove(addedAgainstExisting)
            os.remove(existingAgainstAdded)

        # BLAST added against added.
        addedAgainstAdded = parsedPDB + '/AddedAgainstAdded.txt'
        sequence_BLAST(addedChainsFASTA, addedAgainstAdded, addedDatabaseDir + '/TempDB', blastExecutables + '/psiblast.exe', 2)
        similarities = processPSIoutput.main(addedAgainstAdded)
        for i in similarities:
            query = i[0].split(',')
            hits = i[1].split(',')
            similarity = similarities[i]['Identity']
            matchLength = similarities[i]['MatchLength']
            pairs = [[i,j] for i in query for j in hits]
            similarityTuples = [(i[0], i[0][:-1], i[1], i[1][:-1], similarity, matchLength) for i in pairs]
            c.executemany('INSERT INTO Similarities VALUES (?, ?, ?, ?, ?, ?)', similarityTuples)
            conn.commit()

        # Remove temporary files.
        os.remove(existingChainsFASTA)
        os.remove(addedChainsFASTA)
        os.remove(addedAgainstAdded)
        shutil.rmtree(existingDatabaseDir)
        shutil.rmtree(addedDatabaseDir)

    ##########################
    # Generate output files. #
    ##########################
    # Extract data for file generation.
    c.execute('SELECT * FROM ProteinInformation')
    protInfoResult = c.fetchall()
    c.execute('SELECT * FROM ChainTypes')
    chainTypeResult = c.fetchall()
    c.execute('SELECT * FROM Similarities')
    similarityResult = c.fetchall()

    # Generate PDB entries file.
    allEntries = set([i[1] for i in protInfoResult])
    writeEntries = open(fileAllPDBEntries, 'w')
    for i in allEntries:
        writeEntries.write(i + '\n')
    writeEntries.close()

    # Generate protein information file.
    writeProteinInfo = open(fileProteinInformation, 'w')
    for i in protInfoResult:
        writeProteinInfo.write('\t'.join([str(j) for j in i]) + '\n')
    writeProteinInfo.close()

    # Generate representative file.
    representatives = {}
    for i in protInfoResult:
        chain = i[0]
        sequence = i[-1]
        if sequence in representatives:
            representatives[sequence].append(chain)
        else:
            representatives[sequence] = [chain]
    writeRepresentatives = open(fileRepresentative, 'w')
    for i in representatives:
        chains = representatives[i]
        reprChain = chains[0]
        for j in chains:
            writeRepresentatives.write(j + '\t' + reprChain + '\n')
    writeRepresentatives.close()

    # Generate chain type file.
    writeChains = open(fileChainType, 'w')
    for i in chainTypeResult:
        writeChains.write('\t'.join([str(j) for j in i]) + '\n')
    writeChains.close()

    # Generate similarities file.
    writeSimilarities = open(fileSimilarity, 'w')
    for i in similarityResult:
        writeSimilarities.write('\t'.join([str(j) for j in i]) + '\n')
    writeSimilarities.close()

    # Close the connection to the database.
    c.close()

def sequence_BLAST(inputFile, outputFile, database, BLASTLoc, cores):
    """Will perform the process of BLAST -> PROCESS OUTPUT on inputFile.

    @param inputFile: The FASTA file which needs to be submitted to PSI-BLAST.
    @type inputFile: string
    @param outputFile: The location to write the results of the BLASTing.
    @type outputFile: string
    @param database: The database to BLAST the inputFile protein against
    @param database: string
    @param BLASTLoc: The location of the PSI-BLAST executable.
    @type BLASTLoc: string
    @param cores: The number of threads to create to run BLAST with.
    @type cores: character

    """

    # Setup the parameters for the BLASTing.
    argsPSI = []
    argsPSI.append(BLASTLoc)
    argsPSI.append('-query')
    argsPSI.append(inputFile)
    argsPSI.append('-out')
    argsPSI.append(outputFile)
    argsPSI.append('-evalue')
    argsPSI.append('1')
    argsPSI.append('-num_iterations')
    argsPSI.append('3')
    argsPSI.append('-gap_trigger')
    argsPSI.append('18')
    argsPSI.append('-num_descriptions')
    argsPSI.append('10000')
    argsPSI.append('-num_alignments')
    argsPSI.append('10000')
    argsPSI.append('-dbsize')
    argsPSI.append('0')
    argsPSI.append('-db')
    argsPSI.append(database)
    argsPSI.append('-outfmt')
    argsPSI.append('7 qseqid sseqid pident length evalue')
    argsPSI.append('-num_threads')
    argsPSI.append(str(cores))
    # Perform the BLASTing.
    subprocess.call(argsPSI, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))