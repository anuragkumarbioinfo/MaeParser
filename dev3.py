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
