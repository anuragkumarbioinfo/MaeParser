import sys
mae=sys.argv[1]
atom=sys.argv[2]

from StringIO import StringIO
import csv




class MAEPARSER:
    def __init__(self, mae):
        self.mae = mae
        self.read = ['atom_index','i_m_mmod_type', 'r_m_x_coord', 'r_m_y_coord', 'r_m_z_coord', 'i_m_residue_number', 's_m_mmod_res', 's_m_chain_name', 'i_m_color', 'r_m_charge1', 
                     'r_m_charge2', 's_m_pdb_residue_name', 's_m_pdb_atom_name', 'i_m_atomic_number', 'i_m_formal_charge', 'i_m_representation', 'i_m_visibility', 's_m_color_rgb',
                     's_m_atom_name', 'i_m_secondary_structure', 's_m_label_format', 'i_m_label_color', 's_m_label_user_text', 'r_m_pdb_occupancy', 'i_i_constraint', 'i_i_internal_atom_index',
                     'i_m_pdb_convert_problem', 'i_pdb_PDB_serial', 'i_ppw_het', 's_ppw_CCD_assignment_status', 'r_m_pdb_tfactor', 'i_m_minimize_atom_index', 'i_pa_atomindex',
                     'i_pdb_seqres_index', 's_pa_state', 'i_ppw_water', 'x1', 'x2', 'x3', 'x4', 'x5']
        self.ffx = {'91': 'Molybdenum 6+', '90': 'Molybdenum 5+', '93': 'Lithium neutral', '24': 'Nitrogen - SP', '25': 'N - SP2', '26': 'N - SP3', '27': 'United atom NH  - sp3', '20': 'Oxonium ion - sp2', '21': 'Oxonium ion - sp3', '22': 'EMPTY_TYPE_2', '23': 'Any Oxygen', '95': 'Arsenic Tetrahedral', '28': 'United atom NH2  - sp3', '29': 'United atom NH  - sp2', '94': 'Magnesium neutral', '344': '-CH2-CH(CH3)-CH< : [*][C]([H])([H])[C]([H])([C]([H])([H])[H])[C]([H])([*])[*]', '345': '=C-C(CH3)(CH2-)CH< : [*][C](=[*])[C]([C]([H])([H])[H])([C]([H])([H])[*])[C]([H])([*])[*]', '340': '-CH2-CH2-CH< : [*][C]([H])([H])[C]([H])([H])[C]([H])([*])[*]', '4': 'United atom CH -  sp3', '342': '-CH2-CH2-CH= : [*][C]([H])([H])[C]([H])([H])[C]([H])=[*]', '343': '-CH2-CH2-C(CH3)< : [*][C]([H])([H])[C]([H])([H])[C]([C]([H])([H])[H])([*])[*]', '5': 'United atom CH2 -  sp3', '368': 'NH3-CH2-CH2- : [H][N]([H])([H])[C]([H])([H])[C]([H])([H])[*]', '59': 'Iodine', '58': 'Bromine', '55': 'Boron anion, tetrahedral', '54': 'Boron, trigonal planar', '57': 'Chlorine', '56': 'Fluorine', '51': 'Thiolate anion', '50': 'United atom SH', '53': 'Phosphorus, trivalent', '52': 'Any Sulfur', '371': '(CH3)NH2- : [H][N]([H])([C]([H])([H])[H])[*]', '370': '-CH2-CH(NH3)-CH2- : [*][C]([H])([H])[C]([H])([N]([H])([H])[H])[C]([H])([H])[*]', '373': '(CH3)(CH3)HN- : [H][N]([C]([H])([H])[H])([C]([H])([H])[H])[*]', '372': '(CH3)NH2-CH2-CH2- : [H][N]([H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[*]', '375': '(CH3)(CH3)(CH3)N- : [H][C]([H])([H])[N]([C]([H])([H])[H])([C]([H])([H])[H])[*]', '374': '(CH3)(CH3)HN-CH2-CH2- : [H][N]([C]([H])([H])[H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[*]', '376': '(CH3)(CH3)(CH3)N-CH2-CH2- : [H][C]([H])([H])[N]([C]([H])([H])[H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[*]', '319': '-CH2-CH-CH2- : [*][C]([H])([H])[C]([H])([*])[C]([H])([H])[*]', '318': '-CH2-CH(CH3)-CH2- : [*][C]([H])([H])[C]([H])([C]([H])([H])[H])[C]([H])([H])[*]', '313': '(CH3)(CH3)(CH3)C-CH2- : [H][C]([H])([H])[C]([C]([H])([H])[H])([C]([H])([H])[H])[C]([H])([H])[*]', '312': '(CH3)(CH3)CH-CH2-CH2- : [H][C]([C]([H])([H])[H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[*]', '311': '(CH3)(CH3)(CH3)C- : [H][C]([H])([H])[C]([C]([H])([H])[H])([C]([H])([H])[H])[*]', '310': '(CH3)(CH3)CH- : [H][C]([C]([H])([H])[H])([C]([H])([H])[H])[*]', '317': '-CH(CH3)-CH2- : [*][C]([H])([C]([H])([H])[H])[C]([H])([H])[*]', '316': '-CH2-CH2-CH2-CH2- : [*][C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[*]', '315': '-CH2-CH2-CH2- : [*][C]([H])([H])[C]([H])([H])[C]([H])([H])[*]', '314': '-CH2-CH2- : [*][C]([H])([H])[C]([H])([H])[*]', '115': 'Oxygen anion (-2)', '114': 'Sulfur sulfide anion (-2)', '88': 'Molybdenum 3+', '89': 'Molybdenum 4+', '111': 'Phosphorus cation, tetravalent', '110': 'Sulfur, hexavalent octohedral', '113': 'Sulfur, hexavalent tetrahedral', '112': 'Selenium neutral', '82': 'Cobalt 3+', '83': 'Nickel 2+', '80': 'Iron 3+', '81': 'Cobalt 2+', '86': 'Copper 2+', '87': 'Zinc 2+', '84': 'Nickel 3+', '85': 'Copper +', '108': 'Any Phosphorus', '3': 'Carbon - sp3', '7': 'United atom CH -  sp2', '369': '-CH2-NH2-CH2- : [*][C]([H])([H])[N]([H])([H])[C]([H])([H])[*]', '366': 'NH3- : [H][N]([H])([H])[*]', '367': 'NH3-CH2- : [H][N]([H])([H])[C]([H])([H])[*]', '360': 'CH2(OH)-CH(NH3)-CH(OH)- : [H][C]([O][H])([H])[C]([H])([N]([H])([H])[H])[C]([H])([O][H])[*]', '308': 'CH3-CH2-CH2 : [H][C]([H])([H])[C]([H])([H])[C]([H])([H])[*]', '309': 'CH3-CH2-CH2-CH2- : [H][C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[*]', '301': 'Water: [H][O][H]', '302': 'Hydroniuim ion: [H][O]([H])[H]', '303': 'Hydroxyl ion: [O][H]', '304': 'Place holder for future coarse grained site', '306': 'CH3- : [H][C]([H])([H])[*]', '307': 'CH3-CH2 : [H][C]([H])([H])[C]([H])([H])[*]', '381': '-O-P(=O)(O)-O- : [*][O][P](=[O])([O])[O][*]', '382': '-CH2-O-P(=O)(O)-O- : [*][C]([H])([H])[O][P](=[O])([O])[O][*]', '109': 'Sulfur, tetravalent', '384': '-CH2-C(=O)O- : [*][C]([H])([H])[C](=[O])[O][*]', '385': 'CH(C(=O)O)(NH3)-CH2- : [H][C]([C](=[O])[O])([N]([H])([H])[H])[C]([H])([H])[*]', '102': 'Chloride ion', '103': 'Any Boron', '100': 'Sulfur cation', '101': 'Sulfur sp2 (thioketone)', '106': 'Iodide ion', '107': 'Phosphorus, pentavalent tetrahedral', '104': 'Fluoride ion', '105': 'Bromide ion', '39': 'N- sp2', '38': 'N- sp3', '33': 'United atom NH+ -  sp3', '32': 'N+ - SP3', '31': 'N+ - SP2', '30': 'United atom NH2  - sp2', '37': 'United atom NH2+ - sp2', '36': 'United atom NH+ -  sp2', '35': 'United atom NH3+ - sp3', '34': 'United atom NH2+ - sp3', '383': '-C(=O)O- : [*][C](=[O])[O][*]', '339': '=CH-CH= : [*]=[C]([H])[C]([H])=[*]', '338': '>CH-CH< : [*][C]([H])([*])[C]([H])([*])[*]', '337': 'Phenyl group C6H5 : [*]c1:c([H]):c([H]):c([H]):c([H]):c1([H])', '336': 'Benzene - C6H6 : [H]c1:c([H]):c([H]):c([H]):c([H]):c1([H])', '330': '-C(CH3)=CH- : [*][C]([C]([H])([H])[H])=[C]([H])[*]', '60': 'Silicon', '61': 'Special dummy atom type (FEP)', '62': 'Special Atom Type', '63': 'Lone Pair', '64': 'Any Atom', '65': 'Lithium +', '66': 'Sodium +', '67': 'Potassium +', '68': 'Rubidium +', '69': 'Cesium +', '173': 'Sixteen coordinate', '172': 'Fifteen coordinate', '171': 'Fourteen coordinate', '170': 'Thirteen coordinate', '2': 'Carbon - sp2', '167': 'Ten coordinate', '6': 'United atom CH3 -  sp3', '341': '=CH-CH2-CH< : [*]=[C]([H])[C]([H])([H])[C]([H])([*])[*]', '168': 'Eleven coordinate', '169': 'Icosahedron - twelve coordinate', '8': 'United atom CH2 -  sp2', '164': 'Pentagonal bipyramid - seven coordinate', '165': 'Twisted cube - eight coordinate', '166': 'Nine coordinate', '92': 'Strontium 2+', '160': 'Trigonal bipyramid - five coordinate', '161': 'Octahedral - four coordinate', '162': 'Octahedral - five coordinate', '163': 'Octahedral - six coordinate', '11': 'Carbocation', '10': 'Carbanion', '13': 'EMPTY_TYPE_1', '12': 'Carbon free radical', '15': 'Oxygen - double bond', '14': 'Any Carbon', '17': 'United atom OH', '16': 'Oxygen - single bonds', '19': 'United atom H2O - Water', '18': 'O- (alkoxide, carboxylate) radius picked to match =0 for carboxyls', '391': '-NH-C(=O)- : [*][N]([H])[C](=[O])[*]', '151': 'Isolated atom', '150': 'PI Dummy Atom', '153': 'Linear - two coordinate', '152': 'Linear - single coordinate', '155': 'Trigonal - three coordinate', '154': 'Trigonal - two coordinate', '157': 'Tetrahedral - four coordinate', '156': 'Tetrahedral - three coordinate', '159': 'Trigonal bipyramid - four coordinate', '158': 'Trigonal bipyramid - three coordinate', '48': 'Any Hydrogen', '49': 'Sulfur neutral', '46': 'EMPTY_TYPE_3', '47': 'EMPTY_TYPE_4', '44': 'H-Cation', '45': 'H-Anion', '42': 'H-O(Neut)', '43': 'H-N(Neut)', '40': 'Any Nitrogen', '41': 'H-Electroneutral (e.g. C,S)', '1': 'Carbon - sp', '320': 'General CG site', '326': 'CH2=CH- : [H][C]([H])=[C]([H])[*]', '327': '-CH=CH- : [*][C]([H])=[C]([H])[*]', '9': 'United atom CH -  sp', '328': '-CH=CH-CH2- : [*][C]([H])=[C]([H])[C]([H])([H])[*]', '329': '-CH2-CH=CH-CH2- : [*][C]([H])([H])[C]([H])=[C]([H])[C]([H])([H])[*]', '77': 'Manganese 6+', '76': 'Manganese 5+', '75': 'Manganese 4+', '74': 'Manganese 3+', '73': 'Manganese 2+', '72': 'Magnesium 2+', '71': 'Barium 2+', '70': 'Calcium 2+', '79': 'Iron 2+', '78': 'Manganese 7+', '357': 'CH2(OH)-CH(OH)-CH2- : [H][C]([O][H])([H])[C]([H])([O][H])[C]([H])([H])[*]', '356': 'glycerol - CH2(OH)-CH(OH)-CH2(OH) : [H][C]([O][H])([H])[C]([H])([O][H])[C]([H])([O][H])[H]', '355': '-CH2-CH(OH)-CH2-CH2- : [*][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[C]([H])([H])[*]', '354': '-CH2-CH(OH)-CH2- : [*][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[*]', '353': 'HO-CH2-CH2- : [H][O][C]([H])([H])[C]([H])([H])[*]', '352': 'HO-CH2- : [H][O][C]([H])([H])[*]', '351': 'HO- : [H][O][*]', '359': 'CH2(OH)-CH(OH)-CH(OH)- : [H][C]([O][H])([H])[C]([H])([O][H])[C]([H])([O][H])[*]', '358': 'CH2(OH)-C(OH)-CH2(OH) : [H][C]([O][H])([H])[C]([O][H])([*])[C]([H])([O][H])[H]'}
        self.read2 = {i:n for n,i in enumerate(self.read)}
        self.prop = {}
        self.atomvaluesTAB = {}
        self.atomvalues = {}
        self.res2atm = {}
        self.start()

    def start(self):
        read = False
        with open(self.mae) as filemae:
            counter = 0    
            for lines in filemae:
                if ":::" in lines:
                    counter += 1
                if counter == 5: 
                    read = True
                if counter == 6: 
                    read = False
                if read and ":::" not in lines:
                    if len(lines) > 4:
                        self.load(lines)

    def add(self, key, value):
        name = self.read[key]
        if name not in self.prop:self.prop[name] = []
        self.prop[name].append(value)

    def load(self, line):
        data = StringIO(line)
        reader = csv.reader(data, delimiter=' ')
        line = list(reader)[0][2:]
        #print line
        #print len(line), len(self.read), line
        #raise SystemExit
        for n, i in enumerate(self.read):
            if n >= len(line):
                x = "@"
            else:
                x = line[n]
            #print self.read[n], i
            self.add(n,x)
            

    def getkey(self,key):
        name = self.read[key]
        return self.prop[name]

    def getnkey(self,key, n):
        name = self.read[key]
        return self.prop[name][n+1]

    def getAll(self):
        #print self.read[1], self.prop.keys()
        atoms = len(self.prop[self.read[1]])
        for each in range(atoms):
            #print [len(self.prop[r]) for r in self.read]
            xline = "\t".join([self.prop[r][each] for r in self.read])
            self.atomvaluesTAB[self.prop['atom_index'][each]] = xline
            self.atomvalues[self.prop['atom_index'][each]] = xline.split("\t")
            #print self.prop['atom_index'][each]
        #return self.atomvalues

    def runAll(self):
        if self.atomvalues == {}:
            self.getAll()
            if self.atomvalues == {}:
                raise ValueError("No Values found in MAE")


    def printAtom(self,atomnum):
        self.getAll()
        w = self.atomvalues[atomnum]
        for i in range(len(w)):
            if i == "i_m_mmod_type":
               print self.read[i], "=>" ,w[i], self.ffx[w[i]]
            else:
               print self.read[i], "=>" ,w[i]
                 
    def setres2atm(self):
        self.runAll()
        #load residue to atom number:
        for atom in self.atomvalues:
            name = self.atomvalues[atom]
            resid = name[self.read2['i_m_residue_number']].strip()
            elem  = name[self.read2['s_m_pdb_atom_name']].strip()
            if resid not in self.res2atm: self.res2atm[resid] = {} 
            self.res2atm[resid][elem] = name[self.read2['atom_index']].strip()
    
    def search(self, resid, elem):
        resid = str(resid)
        if self.res2atm == {}: self.setres2atm()
        if resid in self.res2atm:
            if elem in self.res2atm[resid]:
                return self.res2atm[resid][elem] 
            else:
                print "Not Found"
                return self.res2atm[resid]

    def qsite(self):
        self.result = []
        tosearch = [[114,"N"],
                    [113,"C"],
                    [114,"C"],
                    [115,"N"],
                    [144,"N"],
                    [143,"C"],
                    [144,"C"],
                    [145,"N"],
                    [147,"N"],
                    [146,"C"],
                    [147,"C"],
                    [148,"N"],
                    [209,"N"],
                    [208,"C"],
                    [209,"C"],
                    [210,"N"],
                    [243,"N"],
                    [242,"C"],
                    [243,"C"],
                    [244,"N"],
                    [246,"N"],
                    [245,"C"],
                    [246,"C"],
                    [247,"N"]]
        for resid, elem in tosearch:
            self.result.append(self.search(resid,elem))
            
        return "\n".join(["qsitehcap {} {}".format(self.result[i],self.result[i+1]) for i in range(0, len(self.result), 2)]) 


a = MAEPARSER(mae)
a.printAtom(atom)


#print a.qsite()

#print "\nFrozen Atoms\n" + "#" *25 + """
#not ((res.num 114,144,147,209,243,246) OR ((res.ptype " FE ") OR (res.ptype "HOH ") OR (res.ptype "UNK ")) ) 
#""" +  "#" *25 
