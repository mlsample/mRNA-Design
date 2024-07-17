
import numpy as np
import matplotlib.pyplot as plt
import Bio
import os
import copy
import oxDNA_analysis_tools.distance 
import oxDNA_analysis_tools
from oxDNA_analysis_tools.UTILS.RyeReader import describe
import oxpy
from oxDNA_analysis_tools.UTILS.boilerplate import Simulation
import RNA_NN_better as RNA_NN2


def generate_color(ids = [], colors = [], zoom_to = None, script=".custom.color.js"):
    """
        takes a list of ids ,
        and colors as [R,G,B] from 0 to 1, 
        you have to ensure the counts match 
        
        by default zoom on first colored element 

        returns the path to the overlay.js
    """
    id_colors = [[id, color] for id, color in zip(ids, colors)]
    l = len(id_colors)
    if zoom_to is None:
        zoom_to  = ids[0]
    template = f"""
    //default all white 
    let all = [...elements].map(e=>e[1]);
    colorElements(new THREE.Color(1,1,1), all);
    
    let id_colors= {id_colors};

    for(let i = 0; i < {l}; i++){{
        api.selectElementIDs([id_colors[i][0]]);
        colorElements(new THREE.Color(id_colors[i][1][0],
                                      id_colors[i][1][1],
                                      id_colors[i][1][2],
                                    ));
    }}
    api.findElement(elements.get({zoom_to}));
    
    """
    with open(script,"w") as file:
        file.write(template)

    return script


def setup_ssorigami_from_files(topo, hblistfile):
    system,strand = oxDNA_analysis_tools.UTILS.RyeReader.strand_describe(topo)
    seq = []
    for m in strand:
        seq.append(m.btype)
    seq = ''.join(seq)
    ssorigami = ssOrigamiParse(seq)
    ssorigami.load_hblist(hblistfile)
    ssorigami.decode_stretches_from_hb_list()
    ssorigami.assign_nupack_energies(ssorigami._hbstretches,37)
    return ssorigami


def print_uncompatible_regions(origami):
    for i in range(len(origami._hbstretches)):
        leftpairs  = [x[0] for x in origami._hbstretches[i] ]
        rightpairs  = [x[1] for x in origami._hbstretches[i] ][::-1]
        leftseq  = ''.join( [origami._seq[x] for x in leftpairs] )
        rightseq  = ''.join( [origami._seq[x] for x in rightpairs] )
        mismatched = False
        for g in range(len(leftseq)):
            j = len(leftseq) - g -1
            if not are_compatible(leftseq[g],rightseq[j],True):
                mismatched = True
        if mismatched:
            print ('not compatible',leftseq,rightseq, i, origami._hbstretches[i])
    print('')


def get_energy_ratio(non_modified_origami, modified_origami, color_threshold=-0.4, verbose=False):

    diffs = []
    tocolor = []

    for index,originalene in enumerate(non_modified_origami._energies):
        new = modified_origami._energies[index]
        mydiff = (originalene - new)/abs(originalene)
        diffs.append( ( originalene - new)/abs(originalene))
        if mydiff < color_threshold:
            if verbose is True:
                print('Large destabilization: ', len( non_modified_origami._hbstretches[index]), originalene, new)
            for hb in non_modified_origami._hbstretches[index]:
                tocolor.append(hb[0])
                tocolor.append(hb[1])

    return diffs, tocolor


def are_compatible(a, b, include_wobble=False):
    '''
    Returns true if bases a and b can form a base pair
    '''
    if a == 'U':
        a = 'T'
    if b == 'U':
        b = 'T'
    
    if ((a == 'C' and b == 'G') or (a == 'G' and b == 'C') or 
        (a == 'A' and b == 'T') or (a == 'T' and b == 'A')):
            return True

    if include_wobble and ((a == 'G' and b == 'T') or (a == 'T' and b == 'G')):
        return True

    return False


def assign_free_energy_score_DNA(segmentA,segmentB,T):
        pass

def generate_mutual_traps(hbs,stiff=2.0):
    template = "{\n type = mutual_trap\n particle = %d \n ref_particle = %d \n stiff = %f \n r0 = 1.2\n PBC=1\n}\n"
    s = ''
    for hb in hbs:
        s = s + template % (hb[0],hb[1],stiff) 
        s = s + template % (hb[1],hb[0],stiff)
       # s = s +  "{\n type = mutual_trap\n particle = %d \n ref_particle = %d \n stiff = %f \n r0 = 1.2\n PBC=1\n}\n" % (p[1],p[0],stiffness)
    return s


def generate_ffs_ops(hbs,name):
    template = '{ \n order_parameter = bond \n name = %s \n' % (name)
    for i, hb in enumerate(hbs):
        template = template + ' pair%d = %d,%d \n' % (i+1,hb[0],hb[1])
    template = template + '} \n'
    return template


def load_sequence_from_file(seqfile):
    seqs = []
    with open(seqfile) as infile:
            for line in infile:
                seqs.append(line.strip().upper())

    return  ''.join(seqs)


    
class ssOrigamiParse:
    def __init__(self,sequence,isRNA = True):
        '''
        Loads sequence from seqfile; One sequence per line, in 5' to 3' order
        '''
        self._seq =  sequence
        self._isRNA = isRNA
    
    def export_stretches_to_files(self,hblists,prefix,stiffness=0.07):
        allops = ''
        outop = prefix+'op.txt'
        for index,hbs in enumerate( hblists):
            outmutual = prefix+'traps.txt%d' % (index)
                    
            mutual = generate_mutual_traps(hbs,stiffness)
            op = generate_ffs_ops(hbs,'segment%d' % (index))
            allops = allops+op
            with open(outmutual,'w') as out:
                out.write(mutual)
            outcondition = prefix+'condition%d' % (index)
            condition = 'condition1 = {\n  segment%d >= %d \n }\n' % (index,len(hbs))
            with open(outcondition,'w') as out:
                out.write(condition)
        with open(outop,'w') as out:
                out.write(allops) 
    
    def load_pair_distances(self,topofile,conffile,pairlists):
        topology,conf = describe(topofile,conffile)
        self._stretch_distances = []
        alldistances = []
        index = -1
        smallest = 900.
        for i,pairs in enumerate(pairlists):
            leftpairs  = [x[0] for x in pairs ]
            rightpairs = [x[1] for x in pairs ]
            distances = oxDNA_analysis_tools.distance.distance([conf],[topology],[leftpairs],[rightpairs])
            ndist = np.mean( np.ndarray.flatten(np.array(distances)))
            if index == -1 or ndist < smallest:
                smallest = ndist
                index = i
            alldistances.append(ndist )
        
        self._stretch_distances = alldistances
        return index
    
    
    def load_pair_free_energy(self,topofile,conffile,pairlists,T,cutoff):
        topology,conf = describe(topofile,conffile)
        self._stretch_distances = []
        alldistances = []
        allenergies = []
        index = -1
        smallest = 900.
     
        for i,pairs in enumerate(pairlists):
            leftpairs  = [x[0] for x in pairs ]
            rightpairs = [x[1] for x in pairs ]
            distances = oxDNA_analysis_tools.distance.distance([conf],[topology],[leftpairs],[rightpairs])
            leftseq  = ''.join( [self._seq[x] for x in leftpairs] )
            rightseq = ''.join( [self._seq[x] for x in rightpairs[::-1]] )
            #print("Trying index ",i, pairs)
            domain_FE = self.assign_free_energy_score_RNA(leftseq,rightseq,T)  
            ndist = np.mean( np.ndarray.flatten(np.array(distances)))
            if ndist == 0:
                ndist = 0.01
            #distance_penalty = 1.987204259 * T * np.log(ndist)

            #FE = domain_FE - distance_penalty
            if ndist < smallest and domain_FE < cutoff:
                smallest = ndist
                index = i
            alldistances.append(ndist )
            allenergies.append(domain_FE)
        self._stretch_distances = alldistances
        self._free_energies = allenergies
        return index
    
    
    def generate_force_and_op_file_for_stretch(self,hbs,prefix,stiffness = 0.07):
        outop = prefix+'op.txt'
        outmutual = prefix+'trap.txt'
        mutual = generate_mutual_traps(hbs,stiffness)
        op = generate_ffs_ops(hbs,'segment')
        with open(outmutual,'w') as out:
            out.write(mutual)
        outcondition = prefix+'condition.txt' 
        condition = 'action = stop_or \ncondition1 = {\n  segment >= %d \n }\n' % (len(hbs))
        with open(outcondition,'w') as out:
                out.write(condition)
        with open(outop,'w') as out:
                out.write(op) 
    
    def restart_run(self,regions,all_ops_file=None):
        self._formed_regions = []
        self._unformed_regions = copy.deepcopy(regions)
        self._regions_to_form = regions
        if not all_ops_file is None:
            ops = ''
            for index,hbs in enumerate( regions):
                op = generate_ffs_ops(hbs,'segment%d' % (index))
                ops = ops + op
            with open(all_ops_file,'w') as outf:
                outf.write(ops)
        
    def get_hb_regions_state(self,analysis_input_file,opfilename):
        s = Simulation(analysis_input_file)
        s.run()
        s.p.join()

        with open(opfilename) as inf:
            last_line = inf.readlines()[-1]
            vals = [int(x) for x in last_line.strip().split()]
            self._regions_hb_state = vals
            return vals
    
    def check_if_region_formed(self,region_index,region_length,threshold = 0.8):
        if self._regions_hb_state[region_index] / float(region_length) > threshold:
            return True
        else:
            return False
        
    def evaluate_and_iterate_folding_run(self,topofile,conffile,runfolder,prefix='FF'):
        nextone = self.load_pair_distances(topofile,conffile,self._unformed_regions)
        self.generate_force_and_op_file_for_stretch(self._unformed_regions[nextone],runfolder+prefix)
        self._formed_regions.append(self._unformed_regions[nextone])
        self._unformed_regions.remove(self._unformed_regions[nextone])
     
    
    def load_regions_and_iterate_folding_run(self,regions,topofile,conffile,runfolder,prefix,analysis_input_file,opfilename):
        self.get_hb_regions_state(analysis_input_file,opfilename)
        self._formed_regions = []
        self._unformed_regions = []
        for i,region  in enumerate(self._regions_to_form):
            if self.check_if_region_formed(i,len(region)):
                self._formed_regions.append(region)
            else:
                self._unformed_regions.append(region) 
        nextone = self.load_pair_distances(topofile,conffile,self._unformed_regions)
        self.generate_force_and_op_file_for_stretch(self._unformed_regions[nextone],runfolder+prefix)
        #self._formed_regions.append(self._unformed_regions[nextone])
        #self._unformed_regions.remove(self._unformed_regions[nextone])
    
    def load_regions_and_iterate_folding_run_at_T(self,regions,topofile,conffile,runfolder,prefix,analysis_input_file,opfilename,T,cutoff,maxforce=10.,stifness=0.07):
        self.get_hb_regions_state(analysis_input_file,opfilename)
        self._formed_regions = []
        self._unformed_regions = []
        for i,region  in enumerate(self._regions_to_form):
            if self.check_if_region_formed(i,len(region)):
                self._formed_regions.append(region)
            else:
                self._unformed_regions.append(region) 
        
        #nextone = self.load_pair_distances(topofile,conffile,self._unformed_regions)
        
        nextone = self.load_pair_free_energy(topofile,conffile,self._unformed_regions,T,cutoff)
        if nextone == -1:
            print("Cannot find anything, change cutoff")
            return False
        else:
            newdistance = self._stretch_distances[nextone]
            print("Distance is %f and expected force %f " % (newdistance, newdistance * stifness))
            newstifness = stifness
            if newdistance * stifness > maxforce:
                newstifness =   maxforce / newdistance
            print("Setting stiffness", newstifness)
            self.generate_force_and_op_file_for_stretch(self._unformed_regions[nextone],runfolder+prefix,newstifness)
            return True
        #self._formed_regions.append(self._unformed_regions[nextone])
        #self._unformed_regions.remove(self._unformed_regions[nextone])
      
    
    
    
    def load_hblist(self,fname):
         '''
         Loads the list of base pairs from a file that was before created by created by gen_hb_list.sh script
         '''
         self._hblist = []
         with open(fname) as inf:
              for line in inf:
                   vals = line.strip().split()
                   if len(vals) == 2:
                        if int(vals[0]) > int(vals[1]):
                             self._hblist.append([int(vals[1]),int(vals[0])])
                        else:
                             self._hblist.append([int(vals[0]),int(vals[1])])
    
    
    def decode_stretches_from_hb_list(self,cutoff=4):
        hblist = sorted(self._hblist)
        starting_point = hblist[0][0]
        complement_point = hblist[0][1]
        hblist_index = 0
        ending_point = np.max(np.ndarray.flatten(np.array(hblist)))
        all_stretches = []
        new_stretch = []
        while hblist_index < len(hblist)-1:
            new_stretch.append(hblist[hblist_index])
            hblist_index += 1
            if starting_point+1 == hblist[hblist_index][0] and complement_point-1 == hblist[hblist_index][1]:
                starting_point += 1
                complement_point -= 1
            else:
                starting_point =  hblist[hblist_index][0] 
                complement_point =  hblist[hblist_index][1]
                if len(new_stretch) >= cutoff:
                    all_stretches.append(new_stretch)
                new_stretch = []
        
        new_stretch.append(hblist[hblist_index])
        if len(new_stretch) >= cutoff:
            all_stretches.append(new_stretch)
            
        self._hbstretches = all_stretches
        return all_stretches

            
    def compare_hb_list_to_stretches(self,bad_threshold = 4):
        bad = []
        self._found_stretches = []
        for s in self._all_stretches:
            new_stretch = []
            [a_start, a_end, b_start, b_end] = s
            end = b_end

            for a in range(a_start,a_end+1):
                new_stretch.append([a,end])
                end -= 1
                
            self._found_stretches.append(sorted(new_stretch))
        
        all_matches = []
        for index,stretch in enumerate(self._found_stretches):
            #find matching stretch in hbtretches:
            matches = []
            for hindex,hstretch in enumerate(self._hbstretches):
                counter = 0
                for bp in stretch:
                    if bp in hstretch:
                        counter += 1
                if counter > 0:
                    matches.append(counter)
            all_matches.append(matches)
        
        for index,match in enumerate(all_matches):
            wrong = len(self._found_stretches[index])
            for m in match:
                wrong -= m
            if wrong > bad_threshold:
                print ("Stretch length %d, overlaps with %d , %d are wrong" % (len(self._found_stretches[index]),len(match),wrong),self._found_stretches[index] )
                bad.append(self._found_stretches[index])

            else:
                print ("Stretch length %d, overlaps with %d , %d are wrong" % (len(self._found_stretches[index]),len(match),wrong))

                    
                
        return bad
            
    
    def print_cogli_colors(self,hbs,outfile):
        N = len(self._seq)
        jsonfile = '{"RMSF (nm)": ['
        for i in range(N):
            if i not in np.ndarray.flatten(np.array(hbs)):
                jsonfile = jsonfile + '1.0, '
            else:
                 jsonfile = jsonfile + '3.0, '
        
        jsonfile = jsonfile + ']}'
        with open(outfile,'w') as out:
            out.write(jsonfile)
        
    def old_compare_hb_list_to_stretches(self):
        '''
        WRONG DO NOT USE!!
        Defines all hb pairs that exist in the system as detected by stretches
        '''
        self._bad_stretches = []
        self._good_stretches = []
        self._all_stretch_hb = []
        self._missed_hb = []
        for stretch in self._all_stretches:
            [a_start, a_end, b_start, b_end] = stretch
            subhblist = []
            counter = 0
            complement = b_end
            isbad = False
            for i in range(a_start,a_end+1):
                    subhblist.append([i,complement])
                    complement -= 1

            for bp in subhblist:
                 if bp not in self._hblist:
                      counter -= 1
                      isbad = True
                      print('%d %d base pair is extra in stretch ' % (bp[0],bp[1]), stretch)
                      
                 else:
                      counter += 1
            
            
            if counter < 0:
                 self._bad_stretches.append(stretch)
            else:
                 self._good_stretches.append(stretch)
            self._all_stretch_hb = self._all_stretch_hb + subhblist
        
        for hb in self._hblist:
            if hb not in self._all_stretch_hb:
                  self._missed_hb.append(hb)
        
        #filter out the ones which overlap
        return len(self._bad_stretches), len(self._good_stretches), len(self._missed_hb)
    
    
    def indices_overlap(self,stretch,all_stretches):
        maxoverlap = -1
        problematic_ones = []
        start1 = stretch[0] 
        end1 = stretch[1]
        sstart1 = stretch[2]
        eend1   = stretch[3]
        
        for s in all_stretches:
            start2 = s[0]
            end2 = s[1]
            
            sstart2 = s[2]
            eend2   = s[3]
            if max(start1, start2) <= min(end1, end2) or  max(sstart1, sstart2) <= min(eend1, eend2) or  max(sstart1, start2) <= min(eend1, end2) or  max(start1, sstart2) <= min(end1, eend2):
                overlap = end2-start2+1
                if overlap > maxoverlap:
                    maxoverlap = overlap
                
                problematic_ones.append(s)
                
        return maxoverlap, problematic_ones


  
    def indices_overlap_better(self,stretch,all_stretches):
        maxoverlap = -1
        problematic_ones = []
        start1 = stretch[0] 
        end1 = stretch[1]
        sstart1 = stretch[2]
        eend1   = stretch[3]
        
        for s in all_stretches:
            start2 = s[0]
            end2 = s[1]
            
            sstart2 = s[2]
            eend2   = s[3]

            if start1 == start2 and end1 == end2 and sstart1 == sstart2 and eend1 == eend2:
                continue
            elif max(start1, start2) <= min(end1, end2) or  max(sstart1, sstart2) <= min(eend1, eend2) or  max(sstart1, start2) <= min(eend1, end2) or  max(start1, sstart2) <= min(end1, eend2):
                overlap = end2-start2+1
                if overlap > maxoverlap:
                    maxoverlap = overlap
                
                problematic_ones.append(s)
                
        return maxoverlap, problematic_ones

    

    def combine_stretches(self):
        pass
    
    def find_all_stretches(self,min_length=4,wobble_bp=False):
        n = len(self._seq)
        matrix = np.zeros((n,n))
        #maxlen = 0
        res_i = 0
        res_j = 0
        include_wobble = wobble_bp

        self._all_stretches = []
        self._all_energies  = []

        # Populate matrix with base pairs
        for i in range(n):
            for j in range(i,n):
                if are_compatible(self._seq[i], self._seq[j], include_wobble):
                    matrix[i,j] = 1
                    
        # Find longest stretch in matrix
        for i in range(n):
            for j in range(i,n):
                sublen = 0
                k = 0
                while i-k>=0 and k+j<n:
                    if matrix[i-k,j+k] == 1:
                        sublen += 1
                        k += 1
                    else:
                        break
                if sublen >= min_length:
                    #we save this stretch
                    #maxlen = sublen
                    res_i = i
                    res_j = j
                    a_start = res_i - sublen + 1
                    a_end = res_i
                    b_start = res_j
                    b_end = res_j + sublen - 1
                    indices = [a_start, a_end, b_start, b_end]
                    self._all_stretches.append(indices)
        
        #let's resolve the stretches that are overlapping  with each ohter options
        reformed_stretches = copy.deepcopy(self._all_stretches)
        for rs in self._all_stretches:
            mylen = 1+ rs[1] - rs[0]
            if rs not in reformed_stretches:
                continue
            maxcompetition,overlaps = self.indices_overlap(rs,reformed_stretches)
            if maxcompetition > -1:
                if mylen < maxcompetition:
                    reformed_stretches.remove(rs)
                else:
                    for o in overlaps:
                        if o != rs:
                            reformed_stretches.remove(o)
        
        self._all_stretches = reformed_stretches    
        return len(self._all_stretches)
    

    def find_all_stretches_without_conflict(self,min_length=4,tolerance=3,wobble_bp=False):
        #returns all streteches that do not overlap with other stretch by more than tolerance
        n = len(self._seq)
        matrix = np.zeros((n,n))
        #maxlen = 0
        res_i = 0
        res_j = 0
        include_wobble = wobble_bp

        self._all_stretches = []
        self._all_energies  = []

        # Populate matrix with base pairs
        for i in range(n):
            for j in range(i,n):
                if are_compatible(self._seq[i], self._seq[j], include_wobble):
                    matrix[i,j] = 1
                    
        # Find longest stretch in matrix
        for i in range(n):
            for j in range(i,n):
                sublen = 0
                k = 0
                while i-k>=0 and k+j<n:
                    if matrix[i-k,j+k] == 1:
                        sublen += 1
                        k += 1
                    else:
                        break
                if sublen >= min_length:
                    #we save this stretch
                    #maxlen = sublen
                    res_i = i
                    res_j = j
                    a_start = res_i - sublen + 1
                    a_end = res_i
                    b_start = res_j
                    b_end = res_j + sublen - 1
                    indices = [a_start, a_end, b_start, b_end]
                    self._all_stretches.append(indices)
        
        #let's resolve the stretches that are overlapping  with each ohter options
        reformed_stretches = copy.deepcopy(self._all_stretches)
        non_overlapping_stretches = []
        for index, rs in enumerate(self._all_stretches):
            mylen = 1+ rs[1] - rs[0]
            if rs not in reformed_stretches:
                continue
            maxcompetition,overlaps = self.indices_overlap_better(rs,reformed_stretches)
            print("Trying ", index, ' out of ', len(self._all_stretches), rs, maxcompetition, overlaps)
            if maxcompetition > -1:
                if mylen < maxcompetition:
                    print("there is biger stretch, removing", mylen, maxcompetition)
                    reformed_stretches.remove(rs)
                else:
                    for o in overlaps:
                        if o != rs:
                            reformed_stretches.remove(o)
                    print("Candidate for adding",mylen,maxcompetition, tolerance)
                    if maxcompetition <= tolerance:
                        non_overlapping_stretches.append(rs)
                        print("added")
            else:
                non_overlapping_stretches.append(rs)
        
        self._all_stretches = non_overlapping_stretches    

        return len(self._all_stretches)
    

    def assign_energies(self,stretches,T):
        #T is in Kelvin
        self._energies = []
        for stretch in stretches:
            leftpairs  = [x[0] for x in stretch ]
            rightpairs  = [x[1] for x in stretch][::-1]
            subseqA  = ''.join( [self._seq[x] for x in leftpairs] )
            subseqB  = ''.join( [self._seq[x] for x in rightpairs] )
            #if self._isRNA:
            energy = self.assign_free_energy_score_RNA(subseqA,subseqB,T)
            #else:
            #energy = self.assign_free_energy_score_DNA(subseqA,subseqB)    
            self._energies.append(energy)
        return self._energies
    

    def assign_nupack_energies(self,stretches,T):
        #T is in Celsius
        self._energies = []
        for stretch in stretches:
            leftpairs  = [x[0] for x in stretch ]
            rightpairs  = [x[1] for x in stretch][::-1]
            subseqA  = ''.join( [self._seq[x] for x in leftpairs] )
            subseqB  = ''.join( [self._seq[x] for x in rightpairs] )
            #if self._isRNA:
            if ( T > 273.15):
                print("Should be using the temperature in Celsius!")
            energy = RNA_NN2.get_nupack_dg(subseqA,subseqB,T)
            #energy = self.assign_free_energy_score_RNA(subseqA,subseqB,T)
            #else:
            #energy = self.assign_free_energy_score_DNA(subseqA,subseqB)    
            self._energies.append(energy)
        return self._energies



