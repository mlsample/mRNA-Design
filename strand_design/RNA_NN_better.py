import sys
import math
import random
import nupack

GAS_CONSTANT = 1.9858775 * 10**(-3)


#DATAFILEPATH = '/home/petr/studium/oxford/rna_prediction_tools/RNAstructure/data_tables/' 
DATAFILEPATH = '/home/petr/old_homes/old_home/studium/oxford/rna_prediction_tools/RNAstructure/data_tables/'
hstackDG = DATAFILEPATH + 'stack.dat'
hstackDH = DATAFILEPATH + 'stack.dh'

num_to_base = ['A','C','G','U']
base_to_num = {'A' : 0, 'C' : 1, 'G' : 2, 'U' : 3}


def get_nupack_dg(seqA,seqB,temp=37):
	mymodel = nupack.Model(material='rna', ensemble='stacking', celsius=temp)
	sA = nupack.Strand(seqA,name='sA')
	sB = nupack.Strand(seqB,name='sB')
	c = nupack.Complex([sA,sB])
	r = nupack.complex_analysis(complexes=[c],model=mymodel,compute=['pfunc'])
	return r[c].free_energy

	
	

def ReadDSandDH(fileDH,fileDG):
	dgT = 37 + 273.15

	
	finalDG = ReadParams(fileDG)
	finalDH = ReadParams(fileDH)
	finalDS = {}
	for key in finalDH.keys():
		finalDS[key] = 1000.0 * (finalDH[key] - finalDG[key]) / dgT


	return finalDS,finalDH

def are_complementary(baseA, baseB):
	numA = base_to_num[baseA]
	numB = base_to_num[baseB]
	if(numA  + numB  == 3 or numA + numB == 5):
		return True
	else:
		return False



def ReadParams(filename):
	fin = open(filename,'r')
	lines = fin.readlines()
	i = 0
	final_en = {}

	while i < (len(lines)):
		line = lines[i]
		if len(line.split() ) == 4 and line.split() == ['Y','Y','Y','Y']:
			inttypesA = lines[i+5]
			inttypesB = lines[i+6]	
			
			itypesA = inttypesA.split()
			itypesB = inttypesB.split()

			vals = {}
			vals[0] = lines[i+8].split()
			vals[1] = lines[i+9].split()
			vals[2] = lines[i+10].split()
			vals[3] = lines[i+11].split()
		

			for j in range(len(itypesA)):	
					basestrtype = itypesA[j][0] + itypesB[j][0]
					for x in range(4):
						for y in range(4):
							strtype = basestrtype + num_to_base[x] + num_to_base[y]			
							energy = vals[x][j*4+y] 
							if energy != '.' and  are_complementary(num_to_base[x], num_to_base[y]):
								#print strtype, ' = ', energy
								final_en[strtype] = float(energy)

		 
			i = i + 12					
			#sys.exit(1)	
		else:
			i += 1	
	#print final_en
	return final_en			


class NN_calculator:
    def __init__(self):
        self.DS_vals, self.DH_vals = ReadDSandDH(hstackDH,hstackDG)
        
    def get_FEN_DNA(self,seq,compseq,T):
        pass

    def get_FEN_RNA(self,seq,compseq,T):
        length = len(seq)
        terminus = str(seq[0]) + str(compseq[length-1])
        terminusB = str(seq[length-1]) + str(compseq[0])

        DHinitterm = 3.61
        DSinitterm = -1.5

        if( terminus in ['AU','UA','GU','UG'] ) : 
                DHinitterm += 3.72
                DSinitterm += 10.5
                #print >> sys.stderr, ' weak terminal penalty'

        if( terminusB in ['AU','UA','GU','UG'] ) : 
                DHinitterm += 3.72
                DSinitterm += 10.5
                #print >> sys.stderr, ' weak terminal penalty'
        
        DH = DS = 0
        
        for i in range(length-1):
                pair = seq[i] + compseq[length-1-i] + seq[i+1] + compseq[length-1-i-1]
                DH += self.DH_vals[pair]
                DS += self.DS_vals[pair]
                #print >> sys.stderr, 'Adding',pair, DH_vals[pair], DS_vals[pair]
        DH += ( DHinitterm)
        DS += ( DSinitterm) 
        
        return DH*1000 - T*DS

