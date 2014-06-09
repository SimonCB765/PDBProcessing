import sys

import updatelocalPDB
import generateculledsubsets

def main(mmCIFDir, parsedPDB, blastExecutables):
    """Run the updating and culling of the entire PDB.

    :param mmCIFDir:            The directory containing the mmCIF files for the PDB.
    :type mmCIFDir:             string
    :param parsedPDB:           The directory where the results of the parsing and culling will be written.
    :type parsedPDB:            string
    :param blastExecutables:    The location of the BLAST+ executables.
    :type blastExecutables:     string

    """

    updatelocalPDB.main(mmCIFDir, parsedPDB, blastExecutables)
    generateculledsubsets.main(parsedPDB)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])