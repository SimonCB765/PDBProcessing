'''
Created on 11 Jan 2014

@author: Simon Bull
'''

def main(PSIoutput, outputLocation, minAlignLength=20, maxEValue=1.0):
    """Extracts the relevant information from the PSI-BLAST output.

    :param PSIoutput:       The location of the file containing the PSI-BLAST results.
    :type PSIoutput:        string
    :param outputLocation:  The location where the parsed similarity information will be written.
    :type outputLocation:   string
    :param minAlignLength:  The minimum permissible alignment length required for a similarity relationship to be included..
    :type minAlignLength:   int
    :param maxEValue:       The maximum permissible E value required for a similarity relationship to be included.
    :type maxEValue:        float

    """

    currentQuery = ''
    hitsFound = {}
    similaritiesFound = {}

    BLASTOutput = open(PSIoutput, 'r')
    for line in BLASTOutput:
        chunks = line.split()
        if len(chunks) == 0:
            continue
        elif chunks[0] == '#' and chunks[1] == 'Query:':
            # The end of a round has been reached.
            nextQuery = chunks[2]
            if currentQuery != nextQuery:
                # A new query has been found. Write out the hits from the last query.
                for hit in hitsFound:
                    if hit != currentQuery:
                        # Only record the similarity if the query and hit are not the same
                        pair = tuple(sorted([currentQuery, hit]))
                        if not (pair in similaritiesFound and similaritiesFound[pair] >= hitsFound[hit]):
                            # If the pair exists and the recorded similarity is less than the newly found one or the pair does not exist,
                            # then record the new value for the similarity.
                            similaritiesFound[pair] = hitsFound[hit]
            currentQuery = nextQuery
            hitsFound = {}
        elif chunks[0] == currentQuery:
            # An alignment is recorded on the line if the line starts with the query protein.
            hit = chunks[1]
            alignLength = int(chunks[3])
            evalue = float(chunks[4])
            if alignLength >= minAlignLength and evalue <= maxEValue:
                # Only record the hit if the alignment length is long enough and the evalue is large enough.
                hitsFound[hit] = float(chunks[2])
    BLASTOutput.close()

    # Write out the similarity results.
    writeSimilarities = open(outputLocation, 'w')
    writeSimilarities.write('ChainGroupingA\tChainGroupingB\tSimilarity\n')  # Write the header.
    for k, v in similaritiesFound.items():
        writeSimilarities.write(k[0] + '\t' + k[1] + '\t' + str(v) + '\n')
    writeSimilarities.close()