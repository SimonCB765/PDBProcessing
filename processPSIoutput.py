'''
Created on 11 Jan 2014

@author: Simon Bull
'''

def main(PSIoutput, outputLocation):
    """Extracts the relevant information from the PSI-BLAST output.

    @param PSIoutput: The location of the file containing the PSI-BLAST results.
    @type PSIoutput:  string

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
            if alignLength >= 20 and evalue <= 1.0:
                # Only record the hit if the alignment length is long enough and the evalue is large enough.
                hitsFound[hit] = float(chunks[2])
    BLASTOutput.close()

    # Write out the similarity results.
    writeSimilarities = open(outputLocation, 'w')
    writeSimilarities.write('ChainGroupingA\tChainGroupingB\tSimilarity\n')  # Write the header.
    for k, v in similaritiesFound.items():
        writeSimilarities.write(k[0] + '\t' + k[1] + '\t' + str(v) + '\n')
    writeSimilarities.close()


def processWithLimits(PSIoutput):
    """Extracts the relevant information from the PSI-BLAST output.

    @param PSIoutput: The location of the file containing the PSI-BLAST results.
    @type PSIoutput:  string

    """

    currentQuery = ''
    hitsFound = {}
    similaritiesFound = {}

    minPDBID = '00'
    maxPDBID = 'ZZ'
    queryInRange = False
    queryAboveRange = False

    BLASTOutput = open(PSIoutput, 'r')
    for line in BLASTOutput:
        chunks = line.split()
        if len(chunks) == 0:
            continue
        elif chunks[0] == '#' and chunks[1] == 'Query:':
            # If the current line is a comment line with the query protein information on it.
            nextQuery = chunks[2]
            if currentQuery != nextQuery:
                # A new query has been found. Write out the hits from the last query.
                for hit in hitsFound:
                    if hit != currentQuery:
                        # Only record the similarity if the query and hit are not the same
                        pair = tuple(sorted([currentQuery, hit]))
                        if pair in similaritiesFound and similaritiesFound[pair]['Identity'] < hitsFound[hit]['Identity']:
                            # If the pair exists and the recorded similarity is less than the newly found one, then update the record.
                            similaritiesFound[pair]['Identity'] = hitsFound[hit]['Identity']
                            similaritiesFound[pair]['MatchLength'] = hitsFound[hit]['MatchLength']
                        else:
                            # If the pair does not exist, then record the information about it.
                            similaritiesFound[pair] = {'Identity' : hitsFound[hit]['Identity'], 'MatchLength' : hitsFound[hit]['MatchLength']}
            if nextQuery[1:3] >= minPDBID and nextQuery[1:3] <= maxPDBID:
                queryInRange = True
            else:
                queryInRange = False
            if nextQuery[1:3] > maxPDBID:
                queryAboveRange = True
            else:
                queryAboveRange = False
            currentQuery = nextQuery
            hitsFound = {}
        elif chunks[0] == currentQuery:
            # An alignment is recorded on the line if the line starts with the query protein.
            hit = chunks[1]
            if queryInRange and hit[1:3] >= minPDBID:
                # If the query is in the range of interest and the hit is for a pdb ID ordered within or above the range
                hitsFound[hit] = {'Identity' : chunks[2], 'MatchLength' : chunks[3]}
            elif queryAboveRange and (hit[1:3] >= minPDBID and hit[1:3] <= maxPDBID):
                # If the query is above the range and the hit is in the range
                hitsFound[hit] = {'Identity' : chunks[2], 'MatchLength' : chunks[3]}
    BLASTOutput.close()

    return similaritiesFound