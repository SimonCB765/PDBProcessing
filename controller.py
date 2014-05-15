import sys

import updatelocalPDB
import generateculledsubsets

def main(mmCIFDir, parsedPDB, blastExecutables, leafLocation):
	"""
	"""

	updatelocalPDB.main(mmCIFDir, parsedPDB, blastExecutables)
	generateculledsubsets.main(parsedPDB, leafLocation)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])