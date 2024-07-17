import pytest
import strand_design.structring_strand_design as ssd
from pathlib import Path
import shutil
import os
import numpy as np
import ipy_oxdna.dna_structure as dna

class TestStructringStrandDesign:
    """
    Test driven development of the structring strand design algorithm
    """
    def create_dir_with_all_files(self, path):
        file_dir = path / 'test_files'
        file_dir.mkdir()

        ex_structure_file = Path('./strucutres/science_melting/nanobase_files/last_conf_MD2.dat')
        shutil.copy2(ex_structure_file, file_dir)
        
        ex_topology_file = Path('./strucutres/science_melting/nanobase_files/ssRNA_science2.top')
        shutil.copy2(ex_topology_file, file_dir)

        ex_input_md_file = Path('./strucutres/science_melting/nanobase_files/inputMD')
        shutil.copy2(ex_input_md_file, file_dir)

        ex_coding_sequence = Path('./strucutres/science_melting/nanobase_files/coding_seq.txt')
        shutil.copy2(ex_coding_sequence, file_dir)
        
        ex_sequence_dependant_parameters = Path('./strucutres/science_melting/nanobase_files/rna_sequence_dependent_parameters.txt')
        shutil.copy2(ex_sequence_dependant_parameters, file_dir)
        
        return file_dir


    @pytest.fixture
    def dir_all_files(self, tmp_path):
        file_dir = self.create_dir_with_all_files(tmp_path)
        return file_dir


    @pytest.fixture
    def dir_fasta_coding_seq(self, tmp_path):
        file_dir = self.create_dir_with_all_files(tmp_path)
        with open(file_dir / 'coding_seq.txt', 'r') as f:
            lines = f.readlines()
            new_line = '>coding_sequence\n'
            with open(file_dir / 'coding_seq.fasta', 'w') as f:
                f.write(new_line)
                f.writelines(lines)
        os.remove(file_dir / 'coding_seq.txt')

        return file_dir


    @pytest.fixture
    def dir_no_struc(self, tmp_path):
        file_dir = self.create_dir_with_all_files(tmp_path)
        os.remove(file_dir / 'last_conf_MD2.dat')
        return file_dir


    @pytest.fixture
    def dir_no_input_md(self, tmp_path):
        file_dir = self.create_dir_with_all_files(tmp_path)
        os.remove(file_dir / 'inputMD')
        return file_dir


    @pytest.fixture
    def dir_no_coding_sequence(self, tmp_path):
        file_dir = self.create_dir_with_all_files(tmp_path)
        os.remove(file_dir / 'coding_seq.txt')
        return file_dir
    
    
    def set_args(self, my_dir:Path):
        test_args = [
            'structring_strand_design.py',
            '-s', f'{my_dir.as_posix()}/last_conf_MD2.dat',
            '-t', f'{my_dir.as_posix()}/ssRNA_science2.top',
            '-i', f'{my_dir.as_posix()}/inputMD',
            '-c', f'{my_dir.as_posix()}/coding_seq.txt'
            ]
        return test_args
    
    
    def test_collect_args_all_files(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        
        strucutre_file,  topology_file, input_md_file, coding_sequence, output_name = ssd.collect_args()
        assert strucutre_file.exists()
        assert topology_file.exists()
        assert input_md_file.exists()
        assert coding_sequence.exists()
        assert output_name == 'coding_sequence_embbeded_strucutre'
        
        
    def test_collect_args_no_struc(self, dir_no_struc, monkeypatch):
        test_args = self.set_args(dir_no_struc)
        monkeypatch.setattr('sys.argv',test_args)
        
        try:
            ssd.collect_args()
        except ValueError as e:
            assert str(e) == f'Provided strucutre_file path {dir_no_struc}/last_conf_MD2.dat does not exist'
            
            
    def test_collect_args_no_input_md(self, dir_no_input_md, monkeypatch):
        test_args = self.set_args(dir_no_input_md)
        monkeypatch.setattr('sys.argv',test_args)
        
        failed = False
        try:
            ssd.collect_args()
        except ValueError as e:
            assert str(e) == f'Provided input_md_file path {dir_no_input_md}/inputMD does not exist'
            failed = True
        
        assert failed is True  


    def test_collect_args_no_coding_sequence(self, dir_no_coding_sequence, monkeypatch):
        test_args = self.set_args(dir_no_coding_sequence)
        monkeypatch.setattr('sys.argv',test_args)
        
        try:
            strucutre_file,  topology_file, input_md_file, coding_sequence, output_name = ssd.collect_args()
        except ValueError as e:
            assert str(e) == f'Provided coding_sequence path {dir_no_coding_sequence}/coding_seq.txt does not exist'


    def test_parse_coding_sequence_txt(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        assert coding_sequence.suffix == '.txt'
        
        coding_seq_str = ssd.parse_coding_sequence(coding_sequence)
        possible_bases = set('ATCGU')
        assert set(coding_seq_str).issubset(possible_bases)
        assert type(coding_seq_str) == str
        

    def test_parse_coding_sequence_fasta(self, dir_fasta_coding_seq, monkeypatch):
        test_args = self.set_args(dir_fasta_coding_seq)
        test_args = [arg if 'coding_seq.txt' not in arg else f'{dir_fasta_coding_seq}/coding_seq.fasta' for arg in test_args]
        monkeypatch.setattr('sys.argv',test_args)
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        assert coding_sequence.suffix == '.fasta'
        
        coding_seq_str = ssd.parse_coding_sequence(coding_sequence)
        possible_bases = set('ATCGU')
        assert set(coding_seq_str).issubset(possible_bases)
        assert type(coding_seq_str) == str


    def test_get_strucutre(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        
        strucutre = ssd.get_structure(strucutre_file, topology_file)
        assert isinstance(strucutre, dna.DNAStructure)
        
        
    def test_get_hb_list(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        cwd = os.getcwd()
        os.chdir(dir_all_files)
        hb_id_1, hb_id_2 = ssd.get_hblist(strucutre_file, topology_file, input_md_file, traj_file)
        os.chdir(cwd)
        assert len(hb_id_1) == len(hb_id_2)
        
        # I want to be able to make sure that no id in hb_id_1 is in hb_id_2
        assert len(set(hb_id_1).intersection(set(hb_id_2))) == 0


    def test_get_halfs(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        cwd = os.getcwd()
        
        strucutre = ssd.get_structure(strucutre_file, topology_file)
        half, indexes_1, indexes_2 = ssd.get_halfs(strucutre)
        
        assert set(indexes_1).difference(set(indexes_2)) == set(indexes_1)


    def test_make_pair_map(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        cwd = os.getcwd()
        os.chdir(dir_all_files)
        hb_id_1, hb_id_2 = ssd.get_hblist(strucutre_file, topology_file, input_md_file, traj_file)
        os.chdir(cwd)
        pair_map = ssd.make_pair_map(hb_id_1, hb_id_2)
        
        assert np.all(np.isin(np.array(list(pair_map.keys())), np.array(list(pair_map.keys()))))


    def test_get_idx_of_coding_complement(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        
        
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        
        coding_seq_str = ssd.parse_coding_sequence(coding_sequence)
        strucutre = ssd.get_structure(strucutre_file, topology_file)
        
        cwd = os.getcwd()
        os.chdir(dir_all_files)
        hb_id_1, hb_id_2 = ssd.get_hblist(strucutre_file, topology_file, input_md_file, traj_file)
        os.chdir(cwd)
        
        half, indexes_1, indexes_2 = ssd.get_halfs(strucutre)
        pair_map = ssd.make_pair_map(hb_id_1, hb_id_2)
        index_to_seq_map = ssd.get_index_to_seq_map(strucutre) 
                
        coding_indexes, coding_complement = ssd.get_idx_of_coding_complement(indexes_1, indexes_2, pair_map, coding_seq_str)       
        
        for code, comp in zip(coding_indexes, coding_complement):
            code_seq = strucutre.get_base(int(code)).base
            comp_seq = strucutre.get_base(int(comp)).base
            
            assert ssd.get_compseq(code_seq) == comp_seq
            
        assert len(coding_indexes) == len(coding_complement)
            

    def test_find_where_coding_binds_coding(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        
        
        strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file = ssd.collect_args()
        
        coding_seq_str = ssd.parse_coding_sequence(coding_sequence)
        strucutre = ssd.get_structure(strucutre_file, topology_file)
        
        cwd = os.getcwd()
        os.chdir(dir_all_files)
        hb_id_1, hb_id_2 = ssd.get_hblist(strucutre_file, topology_file, input_md_file, traj_file)
        os.chdir(cwd)
        
        half, indexes_1, indexes_2 = ssd.get_halfs(strucutre)
        pair_map = ssd.make_pair_map(hb_id_1, hb_id_2)
        index_to_seq_map = ssd.get_index_to_seq_map(strucutre) 
                
        coding_indexes, coding_complement = ssd.get_idx_of_coding_complement(indexes_1, indexes_2, pair_map, coding_seq_str)       
        
        coding_binds_coding_indexes = ssd.find_where_coding_binds_coding(coding_indexes, coding_complement, pair_map)
        
        assert np.all(np.isin(coding_binds_coding_indexes, coding_indexes))


    def test_mutate_coding()


    def test_main(self, dir_all_files, monkeypatch):
        test_args = self.set_args(dir_all_files)
        monkeypatch.setattr('sys.argv',test_args)
        cwd = os.getcwd()
        os.chdir(dir_all_files)
        out = ssd.main()
        os.chdir(cwd)
        assert out == 0
        
        