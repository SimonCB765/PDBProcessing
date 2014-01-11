import sqlite3

def main(organism=None, minRes=0.0, maxRes=3.0, maxRVal=0.5, minLength=None, maxLength=None, includeNonXray=False, includeCAOnly=False):
    """
    """

    # Generate the SQL query.
    query = 'SELECT chain, sequence FROM ProteinInformation WHERE'
    if organism:
        query += ' organism="' + organism + '" AND'
    if not includeNonXray:
        query += ' experimentType="XRAY" AND'
    if not includeCAOnly:
        query += ' alphaCarbonOnly="0" AND'
    query += ' resolution >= "' + str(minRes) + '" AND'
    query += ' resolution <= "' + str(maxRes) + '" AND'
    query += ' rValueObs <= "' + str(maxRVal) + '"'

    # Connect to the SQLite database
    conn = sqlite3.connect('Leaf.db')
    c = conn.cursor()
    c.execute(query)
    result = c.fetchall()

    # Close the connection to the database.
    c.close()

    # Determine representatives.
    representatives = {}
    for i in result:
        chain = i[0]
        sequence = i[-1]
        if sequence in representatives:
            pass  # Only need to record one representative chain.
        else:
            representatives[sequence] = chain

    # Remove chains that are too short or long.
    if minLength and maxLength:
        for i in representatives:
            seqLen = len(i)
            if seqLen < minLength or seqLen > maxLength:
                del representatives[i]
    elif minLength:
        for i in representatives:
            seqLen = len(i)
            if seqLen < minLength:
                del representatives[i]
    elif maxLength:
        for i in representatives:
            seqLen = len(i)
            if seqLen > maxLength:
                del representatives[i]

    # Return the chains found.
    return sorted([representatives[i] for i in representatives])