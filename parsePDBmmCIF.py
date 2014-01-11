import mmCIFparser

def extract_info(mmCIFFile, tokens='all'):
	"""The parsing basically returns a dictionary with the same structure as the mmCIF file. all the top level block names like _entity and _entry
	have an entry in the returned token dictionary. The keys for each of these top level entries are the sub block names like id and type. Each data
	entry (so for example the data recorded at tokenDict['_entry']['id'] is a list. If there was a loop in the block then the list will have more
	than one entry, else it will have one. If there is more than one entry then each data item within the sub blocks is ordered in the same manner.
	For example, if tokenDict[i][j] == ['a', 'b', 'c'] and tokenDict[i][k] == [1, 2, 3], then one data record in the block had i.j == 'a' and i.k == 1.
	Similarly, i.j == 'b' and i.k == 2, and i.j == 'c' and i.k == 3.
	"""

	tokenDict = mmCIFparser.main(mmCIFFile)

	errorMessage = ''
	if tokens == 'all':
		return tokenDict, errorMessage

	tokensFound = tokenDict.keys()
	subDict = {}
	for i in tokens:
		if '.' in i:
			# If this is True, then the token is something like _entity.id not just _entity.
			chunks = i.split('.')
			mainToken = chunks[0]
			subToken = chunks[1]
			if mainToken in tokensFound:
				if subToken in tokenDict[mainToken].keys():
					if mainToken in subDict:
						subDict[mainToken][subToken] = tokenDict[mainToken][subToken]
					else:
						subDict[mainToken] = {}
						subDict[mainToken][subToken] = tokenDict[mainToken][subToken]
				else:
					errorMessage += 'Main token ' + mainToken + ' found, but sub-token ' + subToken + ' from requested token ' + i + ' not found.\n'
			else:
				errorMessage += 'Main token ' + mainToken + ' from requested token ' + i + ' not found.\n'
		else:
			# If this is True, then the token is something like _entity not _entity.id or _entity.type.
			if i in tokensFound:
				subDict[i] = tokenDict[i]
			else:
				errorMessage += 'Main token ' + i + ' not found.\n'
	return subDict, errorMessage

def main(fileToParse):
	parsedData, errorMessage = extract_info(fileToParse, ['_entry.id', '_entity', '_entity_poly', '_exptl.method', '_atom_site.label_atom_id', '_atom_site.label_entity_id', '_refine', '_reflns', '_struct_ref', '_entity_src_gen', '_entity_src_nat', '_pdbx_entity_src_syn'])
	entryID = parsedData['_entry']['id']
	entityRecords = {}

	# Entity ID and description information.
	for i in range(len(parsedData['_entity']['id'])):
		entityID = parsedData['_entity']['id'][i]
		if 'pdbx_description' in parsedData['_entity']:
			description = parsedData['_entity']['pdbx_description'][i]
		else:
			description = ''
		entityRecords[entityID] = {'description' : description}

	# Entity polymer information.
	if '_entity_poly' in parsedData:
		for i in range(len(parsedData['_entity_poly']['entity_id'])):
			entityID = parsedData['_entity_poly']['entity_id'][i]
			chainIDs = (parsedData['_entity_poly']['pdbx_strand_id'][i].replace(' ', '')).split(',')
			entityRecords[entityID]['chains'] = chainIDs
			polyType = parsedData['_entity_poly']['type'][i]
			if polyType in ['polypeptide(L)', 'polypeptide(D)']:
				entityRecords[entityID]['type'] = 'Protein'
				entityRecords[entityID]['sequence'] = parsedData['_entity_poly']['pdbx_seq_one_letter_code_can'][i]
			else:
				type = []
				if 'polydeoxyribonucleotide' in polyType:
					type.append('DNA')
				if 'polyribonucleotide' in polyType:
					type.append('RNA')
				if 'polysaccharide' in polyType:
					type.append('polysaccharide')
				if type == []:
					entityRecords[entityID]['type'] = 'other'
				else:
					entityRecords[entityID]['type'] = '/'.join(type)

	# Entity source information.
	if '_entity_src_nat' in parsedData:
		for i in range(len(parsedData['_entity_src_nat']['entity_id'])):
			entityID = parsedData['_entity_src_nat']['entity_id'][i]
			scientificName = parsedData['_entity_src_nat']['pdbx_organism_scientific'][i]
			entityRecords[entityID]['scientificName'] = scientificName
	elif '_entity_src_gen' in parsedData:
		for i in range(len(parsedData['_entity_src_gen']['entity_id'])):
			entityID = parsedData['_entity_src_gen']['entity_id'][i]
			scientificName = parsedData['_entity_src_gen']['pdbx_gene_src_scientific_name'][i]
			entityRecords[entityID]['scientificName'] = scientificName
	elif '_pdbx_entity_src_syn' in parsedData:
		for i in range(len(parsedData['_pdbx_entity_src_syn']['entity_id'])):
			entityID = parsedData['_pdbx_entity_src_syn']['entity_id'][i]
			scientificName = parsedData['_pdbx_entity_src_syn']['organism_scientific'][i]
			entityRecords[entityID]['scientificName'] = scientificName

	# External reference (UniProt) information.
	if '_struct_ref' in parsedData:
		for i in range(len(parsedData['_struct_ref']['entity_id'])):
			entityID = parsedData['_struct_ref']['entity_id'][i]
			if 'db_name' in parsedData['_struct_ref']:
				dbName = parsedData['_struct_ref']['db_name'][i]
			else:
				dbName = ''
			if dbName != '':
				entityRecords[entityID]['dbName'] = dbName
				if 'db_code' in parsedData['_struct_ref']:
					dbCode = parsedData['_struct_ref']['db_code'][i]
				else:
					dbCode = ''
				entityRecords[entityID]['dbCode'] = dbCode
			else:
				entityRecords[entityID]['dbName'] = ''
				entityRecords[entityID]['dbCode'] = ''

	# Experimental type.
	expTypes = {'NEUTRON DIFFRACTION' : 'NEUTRON', 'FIBER DIFFRACTION' : 'FIBER', 'X-RAY DIFFRACTION' : 'XRAY', 'ELECTRON MICROSCOPY' : 'EM',
				'FLUORESCENCE TRANSFER' : 'NA', 'POWDER DIFFRACTION' : 'POWDER', 'SOLUTION NMR' : 'NMR', 'INFRARED SPECTROSCOPY' : 'FTIR',
				'ELECTRON CRYSTALLOGRAPHY' : 'NA', 'SOLUTION SCATTERING' : 'NA', 'SOLID-STATE NMR' : 'NMR'}
	if '_exptl' in parsedData:
		if parsedData['_exptl']['method'][0] in expTypes:
			experimentalType = expTypes[parsedData['_exptl']['method'][0]]
		else:
			experimentalType = 'NA'
	else:
		experimentalType = 'NA'

	# Resolution and R Factor Information
	if '_refine' in parsedData:
		if 'ls_d_res_high' in parsedData['_refine']:
			try:
				resolution = float(parsedData['_refine']['ls_d_res_high'][0])
			except:
				resolution = 100
		else:
			resolution = 100
		if 'ls_R_factor_obs' in parsedData['_refine']:
			try:
				rFactorObs = float(parsedData['_refine']['ls_R_factor_obs'][0])
			except:
				rFactorObs = 1
		else:
			rFactorObs = 1
		if 'ls_R_factor_R_free' in parsedData['_refine']:
			try:
				rFactorFree = float(parsedData['_refine']['ls_R_factor_R_free'][0])
			except:
				rFactorFree = 1
		else:
			rFactorFree = 1
	elif '_reflns' in parsedData:
		if 'd_resolution_high' in parsedData['_reflns']:
			try:
				resolution = float(parsedData['_reflns']['d_resolution_high'][0])
			except:
				resolution = 100
		else:
			resolution = 100
		rFactorObs = 1
		rFactorFree = 1
	else:
		resolution = 100
		rFactorObs = 1
		rFactorFree = 1

	# Determine structures with only alpha carbon atoms.
	onlyAlphaCarbon = dict([(i, True) for i in parsedData['_entity']['id']])
	for i in range(len(parsedData['_atom_site']['label_entity_id'])):
		if parsedData['_atom_site']['label_atom_id'][i] != 'CA':
			onlyAlphaCarbon[parsedData['_atom_site']['label_entity_id'][i]] = False
	for i in onlyAlphaCarbon.keys():
		if 'sequence' in entityRecords[i]:
			entityRecords[i]['onlyAlphaCarbon'] = onlyAlphaCarbon[i]

	return entryID, entityRecords, experimentalType, resolution, rFactorObs, rFactorFree