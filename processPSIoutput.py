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
    hitsFound = {}
    similaritiesFound = {}

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
            currentQuery = nextQuery
            hitsFound = {}
        elif chunks[0] == currentQuery:
            # An alignment is recorded on the line if the line starts with the query protein.
            hit = chunks[1]
            hitsFound[hit] = {'Identity' : chunks[2], 'MatchLength' : chunks[3]}
    BLASTOutput.close()

    return similaritiesFound