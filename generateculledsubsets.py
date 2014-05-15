import gzip
import os
import sys

def main(parsedPDB, leafLocation):
    """
    """

    # Import the portions of the Leaf culling that are needed.
    sys.path.append(leafLocation)
    import sparsematrix
    import Leafcull

    # Create the directory to hold the culled subsets.
    subsetsDir = parsedPDB + '/CulledSubsets'
    if not os.path.exists(subsetsDir):
        os.mkdir(subsetsDir)

    # Define the combinations of resolution, r value and sequence identity to use.
    resolutionsOfInterest = [1.6, 1.8, 2.0, 2.2, 2.5, 3.0]
    rValuesOfInterest = [0.25, 1.0]
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
        vertexPairs = [[], []]
        fileSimilarity = parsedPDB + '/Similarity.tsv'
        readSimilarity = open(fileSimilarity, 'r')
        readSimilarity.readlines()  # Strip header.
        for line in readSimilarity:
            chunks = (line.strip()).split('\t')
            source = chunks[0]
            target = chunks[1]
            if source in toCull and target in toCull and float(chunks[3]) >= seqIdentity:
                vertexPairs[0].append(toCull[source])
                vertexPairs[1].append(toCull[target])
        readSimilarity.close()

        # Generate the adjacency list.
        chainsToCull = list(toCull.values())
        indexDict = dict((x, index) for index, x in enumerate(chainsToCull))
        adjacent = sparsematrix.SparseMatrix(len(chainsToCull))
        sourceVertices = [indexDict[x] for x in vertexPairs[0]]
        targetVertices = [indexDict[x] for x in vertexPairs[1]]
        adjacent.addlist(sourceVertices, targetVertices)
        adjacent.addlist(targetVertices, sourceVertices)

        # Perform the culling.
        chainsToRemove, chainsToKeep = Leafcull.main(adjacent, chainsToCull)

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