import os
import shutil
import subprocess
import sys

import parsePDBmmCIF
import processPSIoutput

def main(mmCIFDir, parsedPDB, blastExecutables):
    """
    """

    # Define output files in the TSV format expected by the App Engine bulk uploader (TSV with header).
    if not os.path.exists(parsedPDB):
        os.mkdir(parsedPDB)
    fileChains = parsedPDB + '/Chains.tsv'
    fileSimilarity = parsedPDB + '/Similarity.tsv'
    fileAllFasta = parsedPDB + '/AllChains.fasta'
    fileReprFasta = parsedPDB + '/ReprChains.fasta'

    ##################################################################
    # Go through the mmCIF files and extract the desired information #
    ##################################################################
    writeAllFasta = open(fileAllFasta, 'w')
    uniqueSequences = set([])
    mmCIFDirs = os.listdir(mmCIFDir)
    for i in mmCIFDirs:
        mmCIFFolder = mmCIFDir + '/' + i
        mmCIFFiles = os.listdir(mmCIFFolder)  # Get the files from each subfolder.

        # Process each mmCIF file.
        for j in mmCIFFiles:
            currentFile = mmCIFFolder + '/' + j
            entryID, entityRecords, experimentalType, resolution, rFactorObs, rFactorFree = parsePDBmmCIF.main(currentFile)  # Parse the file.

            # For each record in the entry, examine the data about it. Each entry can have one or more record depending on the different chains
            # recorded in the entry.
            for j in entityRecords.keys():
                if 'type' in entityRecords[j]:
                    # If the record contains type information than examine it further. Only those records with type information are of interest.
                    chains = [entry + chain for entry in entryID for chain in entityRecords[j]['chains']]  # The chains in the record.
                    type = entityRecords[j]['type'].strip()  # The type of the record.

                    if type == 'Protein':
                        # Only interested in the record if it's a protein.
                        dbCode = entityRecords[j]['dbCode'].strip() if 'dbCode' in entityRecords[j] else ''  # External database identifier.
                        if dbCode in ['?', '.']:
                            dbCode = ''
                        dbName = entityRecords[j]['dbName'].strip() if 'dbName' in entityRecords[j] else ''  # External database name.
                        if dbName in ['?', '.']:
                            dbName = ''
                        description = entityRecords[j]['description'].strip() if 'description' in entityRecords[j] else ''  # Record description.
                        if description in ['?', '.']:
                            description = ''
                        onlyAlphaCarbon = entityRecords[j]['onlyAlphaCarbon']  # Whether the structure for the record contains only alpha carbons.
                        scientificName = entityRecords[j]['scientificName'].strip() if 'scientificName' in entityRecords[j] else ''  # Scientific name of the organism the chain belongs to.
                        if scientificName in ['?', '.']:
                            scientificName = ''
                        sequence = entityRecords[j]['sequence'].upper()  # Sequence of the chain.
                        if sequence.count('X') / float(len(sequence)) < 0.5:
                            # If at least 50% of the amino acids in the chain are X, then the 'protein' is deemed to not be a protein.
                            uniqueSequences.add(sequence)
                            for k in chains:
                                # Record the data about each chain.
                                writeAllFasta.write('>' + k + '\t' + str(len(sequence)) + '\t' + experimentalType + '\t' + str(resolution) + '\t' +
                                                    str(rFactorObs) + '\t' + str(rFactorFree) + '\t' + ('no' if onlyAlphaCarbon == 0 else 'yes') + '\t' +
                                                    description + '\t<' + dbName + ' ' + dbCode + '>\t[' + scientificName + ']\n' + sequence + '\n')

    writeAllFasta.close()

    ####################################
    # Determine sequences for BLASTing #
    ####################################
    uniqueSequences = dict((i, index) for index, i in enumerate(uniqueSequences))
    sequencesUsed = set([])
    writeChains = open(fileChains, 'w')
    writeChains.write('\t'.join(['Chain', 'Res', 'RVal', 'SeqLen', 'NonXRay', 'AlphaCarbonOnly', 'ReprGroup']) + '\n')  # Write the header for the chains file.
    writeReprFasta = open(fileReprFasta, 'w')
    readAllFasta = open(fileAllFasta, 'r')
    while True:
        # Read the file two lines at a time.
        identifierLine = readAllFasta.readline().strip()[1:]  # Strip off any whitespace and the > at the front.
        sequence = readAllFasta.readline().strip()
        if not sequence:
            # Reached the end of the file when there is no second line.
            break

        # Write the chain information into the App Engine TSV file.
        sequenceGrouping = uniqueSequences[sequence]
        chunks = identifierLine.split('\t')
        chain = chunks[0]
        resolution = chunks[3]
        rVal = chunks[4]
        nonXRay = 'no' if chunks[2] == 'XRAY' else 'yes'
        alphaCarbonOnly = chunks[6]
        writeChains.write(chain + '\t' + resolution + '\t' + rVal + '\t' + str(len(sequence)) + '\t' + nonXRay + '\t' + alphaCarbonOnly + '\t' +
                          str(sequenceGrouping) + '\n')

        # Record the representative fasta file for BLASTing.
        if not sequenceGrouping in sequencesUsed:
            # Only record the chain if no other chain from its sequence grouping has been recorded.
            writeReprFasta.write('>' + str(sequenceGrouping) + '\n' + sequence + '\n')

        # Update the record of the sequence groupings that have been recorded in the BLASTing file.
        sequencesUsed.add(sequenceGrouping)
    readAllFasta.close()
    writeReprFasta.close()
    writeChains.close()

    ################
    # Run BLASTing #
    ################
    # Generate BLAST database.
    databaseDir = parsedPDB + '/BLASTdatabase'
    if os.path.exists(databaseDir):
        shutil.rmtree(databaseDir)
    os.mkdir(databaseDir)
    makeDBArgs = [blastExecutables + '/makeblastdb', '-in', fileReprFasta, '-out', databaseDir + '/TempDB', '-dbtype', 'prot']
    subprocess.call(makeDBArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # BLAST chains.
    resultsBLAST = parsedPDB + '/ResultsBLAST.txt'
    if os.path.exists(resultsBLAST):
        os.remove(resultsBLAST)
    sequence_BLAST(fileReprFasta, resultsBLAST, databaseDir + '/TempDB', blastExecutables + '/psiblast.exe', 2)
    processPSIoutput.main(resultsBLAST, fileSimilarity)


def sequence_BLAST(inputFile, outputFile, database, BLASTLoc, cores):
    """Will BLAST a given input file.

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
    subprocess.call(argsPSI)#, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)