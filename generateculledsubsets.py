import gzip
import os
import sqlite3
import sys

import selectchains

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
        includeNonXray = False
        includeCAOnly = False

        # Determine the chains that meet the criteria.
        if resolution == 100.0:
            includeNonXray = True
            includeCAOnly = True
        selectedChains = selectchains.main(maxRes=resolution, maxRVal=rValue, minLength=40, includeNonXray=includeNonXray, includeCAOnly=includeCAOnly)

        # Determine similarities between chains.
        chainString = '\',\''.join(selectedChains)
        conn = sqlite3.connect('Leaf.db')
        c = conn.cursor()
        c.execute('SELECT chainA, chainB FROM Similarities WHERE chainA IN (\'' + chainString + '\') AND chainB IN (\'' + chainString + '\') AND similarity >= ?', (seqIdentity,))
        result = c.fetchall()
        c.close()

        # Generate the adjacency list.
        indexDict = dict((selectedChains[x], x) for x in range(len(selectedChains)))
        adjacent = sparsematrix.SparseMatrix(len(selectedChains))
        vertexPairs = list(zip(*result))
        sourceVertices = [indexDict[x] for x in vertexPairs[0]]
        targetVertices = [indexDict[x] for x in vertexPairs[1]]
        adjacent.addlist(sourceVertices, targetVertices)
        adjacent.addlist(targetVertices, sourceVertices)

        # Perform the culling.
        chainsToCull, chainsToKeep = Leafcull.main(adjacent, selectedChains)

        # Get data on the kept chains.
        conn = sqlite3.connect('Leaf.db')
        c = conn.cursor()
        c.execute('SELECT * FROM ProteinInformation WHERE chain IN (\'' + '\',\''.join(chainsToKeep) + '\')')
        keptChainInfo = c.fetchall()
        c.close()

        # Write out the results.
        xrayCAInfo = ''
        if includeNonXray:
            xrayCAInfo = '_INCLNONXRAY'
        if includeCAOnly:
            xrayCAInfo += '_INCLCAONLY'
        outputLocation = subsetsDir + '/SeqIden_' + str(seqIdentity) + '_Res_' + str(resolution) + '_RVal_' + str(rValue) + xrayCAInfo + '.fasta.gz'
        with gzip.open(outputLocation, 'w') as writeKept:
            for i in keptChainInfo:
                sequence = i[-1]
                writeKept.write(bytes('>' + i[0] + '\t' + str(len(sequence)) + '\t' + i[2] + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i[5]) + '\t' +
                                      ('no' if i[6] == 0 else 'yes') + '\t' + i[7] + '\t<' + i[8] + ' ' + i[9] + '>\t[' + i[10] + ']\n' + sequence + '\n',
                                      'UTF-8'))