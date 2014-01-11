'''
Created on 11 Jan 2014

@author: Simon Bull
'''

def main(PSIoutput):
    """Extracts the relevant information from the PSI-BLAST output.

    @param PSIoutput: The location of the file containing the PSI-BLAST results.
    @type PSIoutput:  string

    """

    currentQuery = ''
    driftFound = False
    hitsFoundLastIteration = {}
    hitsFoundThisIteration = {}
    similaritiesFound = {}

    BLASTOutput = open(PSIoutput, 'r')
    for line in BLASTOutput:
        chunks = line.split()
        if len(chunks) == 0:
            # If the line is a blank line, then the end of a query has been reached.
            driftFound = False
            for hit in hitsFoundThisIteration:
                if hit != currentQuery:
                    # Only record the similarity if the query and hit are not the same
                    pair = tuple(sorted([currentQuery, hit]))
                    if pair in similaritiesFound and similaritiesFound[pair]['Identity'] < hitsFoundThisIteration[hit]['Identity']:
                        # If the pair exists and the recorded similarity is less than the newly found one, then update the record.
                        similaritiesFound[pair]['Identity'] = hitsFoundThisIteration[hit]['Identity']
                        similaritiesFound[pair]['MatchLength'] = hitsFoundThisIteration[hit]['MatchLength']
                    else:
                        # If the pair does not exist, then record the information about it.
                        similaritiesFound[pair] = {'Identity' : hitsFoundThisIteration[hit]['Identity'], 'MatchLength' : hitsFoundThisIteration[hit]['MatchLength']}
            hitsFoundLastIteration = {}
            hitsFoundThisIteration = {}
        elif driftFound:
            continue
        elif chunks[0] == '#' and chunks[1] == 'Query:':
            # If the current line is a comment line with the query protein information on it.
            for key in hitsFoundLastIteration:
                # Compare currentRoundHits with previousRoundHits to check for drift.
                if not key in hitsFoundThisIteration:
                    # If a hit is recorded in hitsFoundLastIteration but not in hitsFoundThisIteration then drift is deemed to have occurred.
                    hitsFoundThisIteration = hitsFoundLastIteration
                    driftFound = True
            if not driftFound:
                currentQuery = chunks[2]
                hitsFoundLastIteration = hitsFoundThisIteration
                hitsFoundThisIteration = {}
        elif chunks[0] == currentQuery:
            # An alignment is recorded on the line if the line starts with the query protein.
            hit = chunks[1]
            hitsFoundThisIteration[hit] = {'Identity' : chunks[2], 'MatchLength' : chunks[3]}
    BLASTOutput.close()

    return similaritiesFound