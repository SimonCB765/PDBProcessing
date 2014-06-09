import gzip
import os
import Leafcull

def main(parsedPDB):
    """Cull the entire PDB at different quality criterion.

    :param parsedPDB:   The directory where the results of the culling will be written.
    :type parsedPDB:    string

    """

    # Create the directory to hold the culled subsets.
    subsetsDir = parsedPDB + '/CulledSubsets'
    if not os.path.exists(subsetsDir):
        os.mkdir(subsetsDir)

    # Define the combinations of resolution, r value and sequence identity to use.
    resolutionsOfInterest = [1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0]
    rValuesOfInterest = [0.25, 0.5, 0.75, 1.0]
    sequenceIdentities = [20, 25, 30, 40, 50, 60, 70, 80, 90]
    tuplesToDo = [(i, j, k) for i in resolutionsOfInterest for j in rValuesOfInterest for k in sequenceIdentities]
    tuplesToDo.extend([(100.0, 1.0, i) for i in sequenceIdentities])
    tuplesToDo = set(tuplesToDo)

    # Generate the culled lists.
    for i in tuplesToDo:
        resolution = i[0]
        rValue = i[1]
        seqIdentity = i[2]
        includeNonXrayAndCAOnly = 1 if resolution == 100.0 else 0

        # Determine the chains that meet the criteria.
        toCull = {}
        fileChains = parsedPDB + '/Chains.tsv'
        readChains = open(fileChains, 'r')
        readChains.readline()  # Strip the header.
        for line in readChains:
            chunks = (line.strip()).split('\t')
            chain = chunks[0]
            res = float(chunks[1])
            rVal = float(chunks[2])
            seqLen = int(chunks[3])
            nonXRay = 1 if chunks[4] == 'yes' else 0
            alphaCarbonOnly = 1 if chunks[5] == 'yes' else 0
            reprGroup = chunks[6]
            if (res <= resolution) and (rVal <= rValue) and (seqLen >= 40) and (nonXRay <= includeNonXrayAndCAOnly) and (alphaCarbonOnly <= includeNonXrayAndCAOnly):
                toCull[reprGroup] = chain
        readChains.close()

        # Determine similarities between representative groups that need culling.
        adjList = {}
        fileSimilarity = parsedPDB + '/Similarity.tsv'
        readSimilarity = open(fileSimilarity, 'r')
        readSimilarity.readlines()  # Strip header.
        for line in readSimilarity:
            chunks = (line.strip()).split('\t')
            chainA = chunks[0]
            chainB = chunks[1]
            if chainA in toCull and chainB in toCull and float(chunks[3]) >= seqIdentity:
                # The sequences are in the set to be culled and are too similar.
                if chainA in adjList:
                    adjList[chainA].add(chainB)
                else:
                    adjList[chainA] = set([chainB])
                if chainB in adjList:
                    adjList[chainB].add(chainA)
                else:
                    adjList[chainB] = set([chainA])
        readSimilarity.close()

        # Perform the culling.
        chainsToRemove = Leafcull.main(adjList)
        chainsToKeep = [toCull[i] for i in toCull if not i in chainsToRemove]

        # Write out the kept chains.
        xrayCAInfo = '_INCLNONXRAY_INCLCAONLY' if includeNonXrayAndCAOnly else ''
        outputLocation = subsetsDir + '/SeqIden_' + str(seqIdentity) + '_Res_' + str(resolution) + '_RVal_' + str(rValue) + xrayCAInfo + '.fasta.gz'
        with gzip.open(outputLocation, 'w') as writeKept:
            fileAllFasta = parsedPDB + '/AllChains.fasta'
            readAllFasta = open(fileAllFasta, 'r')
            while True:
                # Read the file two lines at a time.
                identifierLine = readAllFasta.readline()
                sequence = readAllFasta.readline()
                if not sequence:
                    # Reached the end of the file when there is no second line.
                    break

                chain = identifierLine[1:6]  # Get the chain identifier.
                if chain in chainsToKeep:
                    writeKept.write(bytes(identifierLine + sequence, 'UTF-8'))
            readAllFasta.close()