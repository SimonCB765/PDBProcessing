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
        seqLen = len(sequence)
        lengthOK = True

        # Determine if the chain meets the length criteria.
        if minLength and maxLength:
            if seqLen < minLength or seqLen > maxLength:
                lengthOK = False
        elif minLength:
            if seqLen < minLength:
                lengthOK = False
        elif maxLength:
            if seqLen > maxLength:
                lengthOK = False

        if lengthOK:
            if sequence in representatives:
                pass  # Only need to record one representative chain.
            else:
                representatives[sequence] = chain

    # Return the chains found.
    return sorted([representatives[i] for i in representatives])