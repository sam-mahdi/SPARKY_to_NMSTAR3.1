import re

"""
This script has various functions that look at each file and check various formatting/labeling errors. I.E.:
If a non-existing amino acid was added (e.g. Z)
If the amino acid or atom as forgotten (e.g. 334N or Y334)
If the atoms were incorrectly typed (N was forgotten or typed incorrectly)
If the wrong name was used for the amino acid (e.g. HA instead of HA2 for Glycine)
If the 2 dimensions don't match (e.g. Y112N-something-E148HN)
If the i-1 is not the i-1 (e.g. Y112N-F334N-Y112HN)
If there is something else wrong, or if the format is so entirely off that the script cannot continue, an error will pop up to inform the user which peak should be corrected
"""

accepted_letters=['A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','V','W','Y']

def NHSQC_checker(nhsqc_file):
    print('Checking NHSQC')
    with open(nhsqc_file) as nhsqc:
      for nhsqc_lines in nhsqc:
          if nhsqc_lines.strip().split() == []:
              continue
          if nhsqc_lines.strip().split()[0] in {'Assignment','?-?'}:
              continue
          nhsqc_split=nhsqc_lines.strip().split()
          try:
              amino_acid=nhsqc_lines.strip().split()[0][0]
              atom_one=(re.search('[A-Z]\d+(\w+)',nhsqc_split[0].split('-')[0])).group(1)
              atom_two=nhsqc_split[0].split('-')[1]
              if amino_acid not in accepted_letters:
                  print(f'Amino Acid {nhsqc_split[0]} amino acid is improperly labeled')
              if re.search('[A-Z]\d+\w+',nhsqc_lines.strip().split()[0]) is None:
                  print(f'Amino Acid {nhsqc_split[0]} format is wrong')
              if atom_one not in {'N','NE2','ND2'}:
                  print(f'Amino Acid {nhsqc_split[0]} amide nitrogen is improperly labeled (NH2 should be NE2/ND2)')
              if atom_two not in {'HN','HE21','HE22','HD21','HD22'}:
                  print(f'Amino Acid {nhsqc_split[0]} amide hydrogen is improperly labeled (NH2 should be labeled HE/D 21/22)')
          except:
              print(f'Program could not analyze, stopped at \n {nhsqc_lines}please check peak, correct, and rerun\n')
    print('\n')

def HNCA_checker(hnca_file):
    print('Checking HCNA')
    with open(hnca_file) as hnca:
        for hnca_lines in hnca:
            if hnca_lines.strip().split() == []:
                continue
            if hnca_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            hnca_split=hnca_lines.strip().split()
            try:
                amino_acid=hnca_lines.strip().split()[0][0]
                if amino_acid not in accepted_letters:
                    print(f'Amino Acid {hnca_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',hnca_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {hnca_split[0]} format is wrong')
                if hnca_split[0].split('-')[1] != 'CA' and hnca_split[0].split('-')[2] == 'HN':
                    print(f"Amino Acid {hnca_split[0]} CA is improperly labeled")
                if hnca_split[0].split('-')[2] != 'HN':
                    if (re.search('^\w+\d+',hnca_split[0].split('-')[0])).group(0) != (re.search('^\w+\d+',hnca_split[0].split('-')[2])).group(0):
                        print(f"Amino Acid {hnca_split[0]} nitrogen or amide is improperly labeled")
                        continue
                    i_atom=re.search('(\d+)(\w+)',hnca_split[0].split('-')[0])
                    i_minus_atom=re.search('(\d+)(\w+)',hnca_split[0].split('-')[1])
                    if int(i_atom.group(1)) != int(i_minus_atom.group(1))+1:
                        print(f"Amino Acid {hnca_split[0]} i-1 is improperly labeled (it is not labeled as the i-1)")
                    if i_minus_atom.group(2) != 'CA':
                        print(f"Amino Acid {hnca_split[0]} CA is improperly labeled (its CA is not labeled properly)")
            except:
                print(f'Program could not analyze, stopped at \n {hnca_lines}please check peak, correct, and rerun\n')
    print('\n')
def HNCACB_checker(hncacb_file):
    print('Checking HNCACB')
    with open(hncacb_file) as hncacb:
        for hncacb_lines in hncacb:
            if hncacb_lines.strip().split() == []:
                continue
            if hncacb_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            hncacb_split=hncacb_lines.strip().split()
            try:
                amino_acid=hncacb_lines.strip().split()[0][0]
                if amino_acid not in accepted_letters:
                    print(f'Amino Acid {hncacb_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',hncacb_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {hncacb_split[0]} format is wrong')
                if hncacb_split[0].split('-')[1] not in {'CA','CB'} and hncacb_split[0].split('-')[2] == 'HN':
                    print(f"Amino Acid {hncacb_split[0]} CA is improperly labeled")
                if hncacb_split[0].split('-')[2] != 'HN':
                    if hncacb_split[0].split('-')[2] == 'N':
                        print(f"Amino Acid {hncacb_split[0]} amide is improperly labeled")
                        continue
                    if (re.search('^\w+\d+',hncacb_split[0].split('-')[0])).group(0) != (re.search('^\w+\d+',hncacb_split[0].split('-')[2])).group(0):
                        print(f"Amino Acid {hncacb_split[0]} nitrogen or amide is improperly labeled")
                        continue
                    i_atom=re.search('(\d+)(\w+)',hncacb_split[0].split('-')[0])
                    i_minus_atom=re.search('(\d+)(\w+)',hncacb_split[0].split('-')[1])
                    if int(i_atom.group(1)) != int(i_minus_atom.group(1))+1:
                        print(f"Amino Acid {hncacb_split[0]} i-1 is improperly labeled (it is not labeled as the i-1)")
                    if i_minus_atom.group(2) not in {'CA','CB'}:
                        print(f"Amino Acid {hncacb_split[0]} CA or CB is not labeled properly")
            except:
                print(f'Program could not analyze, stopped at \n {hncacb_lines}please check peak, correct, and rerun\n')
    print('\n')
def HNCOCA_checker(hncoca_file):
    print('Checking HNCOCA')
    with open(hncoca_file) as hnca1:
        for hnca1_lines in hnca1:
            if hnca1_lines.strip().split() == []:
                continue
            if hnca1_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            hnca1_split=hnca1_lines.strip().split()
            try:
                amino_acid=hnca1_lines.strip().split()[0][0]
                if amino_acid not in accepted_letters:
                    print(f'Amino Acid {hnca1_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',hnca1_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {hnca1_split[0]} format is wrong')
                if hnca1_split[0].split('-')[2] in {'N','HN'}:
                    print(f"Amino Acid {hnca1_split[0]} amide is improperly labeled")
                    continue
                if (re.search('^\w+\d+',hnca1_split[0].split('-')[0])).group(0) != (re.search('^\w+\d+',hnca1_split[0].split('-')[2])).group(0):
                    print(f"Amino Acid {hnca1_split[0]} nitrogen or amide is improperly labeled")
                i_atom=re.search('(\d+)(\w+)',hnca1_split[0].split('-')[0])
                i_minus_atom=re.search('(\d+)(\w+)',hnca1_split[0].split('-')[1])
                if int(i_atom.group(1)) != int(i_minus_atom.group(1))+1:
                    print(f"Amino Acid {hnca1_split[0]} i-1 is improperly labeled (it is not labeled as the i-1)")
                if i_minus_atom.group(2) != 'CA':
                    print(f"Amino Acid {hnca1_split[0]} CA is improperly labeled")
            except:
                print(f'Program could not analyze, stopped at \n {hnca1_lines}please check peak, correct, and rerun\n')
    print('\n')
def HNCO_checker(hnco_file):
    print('Checking HNCO')
    with open(hnco_file) as hnco:
        for hnco_lines in hnco:
            if hnco_lines.strip().split() == []:
                continue
            if hnco_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            hnco_split=hnco_lines.strip().split()
            try:
                amino_acid=hnco_lines.strip().split()[0][0]
                if amino_acid not in accepted_letters:
                    print(f'Amino Acid {hnco_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',hnco_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {hnco_split[0]} format is wrong')
                if hnco_split[0].split('-')[2] in {'N','HN'}:
                    print(f"Amino Acid {hnco_split[0]} amide is improperly labeled")
                    continue
                if (re.search('^\w+\d+',hnco_split[0].split('-')[0])).group(0) != (re.search('^\w+\d+',hnco_split[0].split('-')[2])).group(0):
                    print(f"Amino Acid {hnco_split[0]} nitrogen or amide is improperly labeled")
                i_atom=re.search('(\d+)(\w+)',hnco_split[0].split('-')[0])
                i_minus_atom=re.search('(\d+)(\w+)',hnco_split[0].split('-')[1])
                if int(i_atom.group(1)) != int(i_minus_atom.group(1))+1:
                    print(f"Amino Acid {hnco_split[0]} i-1 is improperly labeled (it is not labeled as the i-1)")
                if i_minus_atom.group(2) != 'C':
                    print(f"Amino Acid {hnco_split[0]} C is improperly labeled")
            except:
                print(f'Program could not analyze, stopped at \n {hnco_lines}please check peak, correct, and rerun\n')
    print('\n')
def CHSQC_checker(chsqc_file):
    print('Checking CHSQC')
    with open(chsqc_file) as chsqc:
        for chsqc_lines in chsqc:
            if chsqc_lines.strip().split() == []:
                continue
            if chsqc_lines.strip().split()[0] in {'Assignment','?-?'}:
                continue
            chsqc_split=chsqc_lines.strip().split()
            try:
                amino_acids=chsqc_lines.strip().split()[0][0]
                if amino_acids not in accepted_letters:
                    print(f'Amino Acid {chsqc_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',chsqc_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {chsqc_split[0]} format is wrong')
                search=re.search('([A-Z])\d+(\w+)-(\w+\d*)',chsqc_split[0])
                if search.group(0) is None:
                    print(f'Amino Acid {chsqc_split[0]} is improperly formatted')
                amino_acid=search.group(1)
                atom_1=search.group(2)
                atom_2=search.group(3)
                if amino_acid == 'A':
                    if atom_1 not in {'CA','CB'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'R':
                    if atom_1 not in {'CA','CB','CG','CD'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid in {'N','D','W','S','F','H','C'}:
                    if atom_1 not in {'CA','CB'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid in {'Q','E'}:
                    if atom_1 not in {'CA','CB','CG'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3','HG2','HG3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'G':
                    if atom_1 not in {'CA'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA2','HA3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'I':
                    if atom_1 not in {'CA','CB','CG1','CG2','CD1'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB','HG12','HG13','HD1','HG2'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'L':
                    if atom_1 not in {'CA','CB','CG','CD1','CD2'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3','HG','HD1','HD2'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'K':
                    if atom_1 not in {'CA','CB','CG','CD','CE'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3','HE2','HE3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'M':
                    if atom_1 not in {'CA','CB','CG','CE'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'P':
                    if atom_1 not in {'CA','CB','CG','CD'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB2','HB3','HG2','HG3','HD2','HD3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'T':
                    if atom_1 not in {'CA','CB','CG2'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB','HG2'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if amino_acid == 'V':
                    if atom_1 not in {'CA','CB','CG1','CG2'}:
                        print(f'Amino Acid {chsqc_split[0]} carbon improperly labeled')
                    if atom_2 not in {'HA','HB','HG1','HG2'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
            except:
                print(f'Program could not analyze, stopped at \n {chsqc_lines}please check peak, correct, and rerun\n')
    print('\n')


def HBHACONH_checker(hbhaconh_file):
    print('Checking HBHACONH')
    with open(hbhaconh_file) as hbhaconh:
        for hbhaconh_lines in hbhaconh:
            if hbhaconh_lines.strip().split() == []:
                continue
            if hbhaconh_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            hbhaconh_split=hbhaconh_lines.strip().split()
            try:
                amino_acids=hbhaconh_lines.strip().split()[0][0]
                if amino_acids not in accepted_letters:
                    print(f'Amino Acid {hbhaconh_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',hbhaconh_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {hbhaconh_split[0]} format is wrong')
                if hbhaconh_split[0].split('-')[2] in {'N','HN'}:
                    print(f"Amino Acid {hbhaconh_split[0]} amide is improperly labeled")
                    continue
                if (re.search('^\w+\d+',hbhaconh_split[0].split('-')[0])).group(0) != (re.search('^\w+\d+',hbhaconh_split[0].split('-')[2])).group(0):
                    print(f"Amino Acid {hbhaconh_split[0]} nitrogen or amide is improperly labeled")
                i_atom=re.search('(\d+)(\w+)',hbhaconh_split[0].split('-')[0])
                i_minus_atom=re.search('([A-Z])(\d+)(\w+)',hbhaconh_split[0].split('-')[1])
                if int(i_atom.group(1)) != int(i_minus_atom.group(2))+1:
                    print(f"Amino Acid {hbhaconh_split[0]} i-1 is improperly labeled (it is not labeled as the i-1)")
                if i_minus_atom.group(1) == 'A':
                    if i_minus_atom.group(3) not in {'HA','HB'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'R':
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) in {'N','D','W','S','F','H','C'}:
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) in {'Q','E'}:
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3','HG3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'G':
                    if i_minus_atom.group(3) not in {'HA2','HA3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'I':
                    if i_minus_atom.group(3) not in {'HA','HB'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'L':
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3','HG'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'K':
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'M':
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'P':
                    if i_minus_atom.group(3) not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'T':
                    if i_minus_atom.group(3) not in {'HA','HB'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
                if i_minus_atom.group(1) == 'V':
                    if i_minus_atom.group(3) not in {'HA','HB'}:
                        print(f'Amino Acid {chsqc_split[0]} hydrogen improperly labeled')
            except:
                print(f'Program could not analyze, stopped at \n {hbhaconh_lines}please check peak, correct, and rerun\n')
    print('\n')
def CCH_TOCSY_checker(cch_tocsy_file):
    print('Checking CCH_TOCSY')
    with open(cch_tocsy_file) as ccc_tocsy:
        for ccc_tocsy_lines in ccc_tocsy:
            if ccc_tocsy_lines.strip().split() == []:
                continue
            if ccc_tocsy_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            ccc_tocsy_split=ccc_tocsy_lines.strip().split()
            try:
                amino_acids=ccc_tocsy_lines.strip().split()[0][0]
                if amino_acids not in accepted_letters:
                    print(f'Amino Acid {ccc_tocsy_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',ccc_tocsy_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {ccc_tocsy_split[0]} format is wrong')
                search=re.search('([A-Z])\d+(\w+)',ccc_tocsy_split[0].split('-')[0])
                amino_acid=search.group(1)
                atom=search.group(2)
                if amino_acid == 'A':
                    if atom not in {'CA','CB'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'R':
                    if atom not in {'CA','CB','CG','CD'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG','CD'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid in {'N','D','W','S','F','H','C'}:
                    if atom not in {'CA','CB'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid in {'Q','E'}:
                    if atom not in {'CA','CB','CG'}  or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG2','HG3'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'G':
                    if atom not in {'CA'} or ccc_tocsy_split[0].split('-')[2] not in {'HA2','HA3'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'I':
                    if atom not in {'CA','CB','CG1','CG2','CD1'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB','HG12','HG13','HD1','HG2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG1','CG2','CD1'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'L':
                    if atom not in {'CA','CB','CG','CD1','CD2'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG','HD1','HD2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG','CD1','CD2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'K':
                    if atom not in {'CA','CB','CG','CD','CE'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3','HE2','HE3'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG','CD','CE'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'M':
                    if atom not in {'CA','CB','CG','CE'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG','CE'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'P':
                    if atom not in {'CA','CB','CG','CD'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG','CD'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'T':
                    if atom not in {'CA','CB','CG2'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB','HG2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
                if amino_acid == 'V':
                    if atom not in {'CA','CB','CG1','CG2'} or ccc_tocsy_split[0].split('-')[2] not in {'HA','HB','HG1','HG2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if ccc_tocsy_split[0].split('-')[1] not in {'CA','CB','CG1','CG2'}:
                        print(f'Amino Acid {ccc_tocsy_split[0]} TOCSY carbon improperly labeled')
            except:
                print(f'Program could not analyze, stopped at \n {ccc_tocsy_lines}please check peak, correct, and rerun\n')
    print('\n')
def HCCH_TOCSY_checker(hcch_tocsy_file):
    print('Checking HCCH_TOCSY')
    with open(hcch_tocsy_file) as hcch_tocsy:
        for hcch_tocsy_lines in hcch_tocsy:
            if hcch_tocsy_lines.strip().split() == []:
                continue
            if hcch_tocsy_lines.strip().split()[0] in {'Assignment','?-?-?'}:
                continue
            hcch_tocsy_split=hcch_tocsy_lines.strip().split()
            try:
                amino_acids=hcch_tocsy_lines.strip().split()[0][0]
                if amino_acids not in accepted_letters:
                    print(f'Amino Acid {hcch_tocsy_split[0]} amino acid is improperly labeled')
                if re.search('[A-Z]\d+\w+',hcch_tocsy_lines.strip().split()[0]) is None:
                    print(f'Amino Acid {hcch_tocsy_split[0]} format is wrong')
                search=re.search('([A-Z])\d+(\w+)',hcch_tocsy_split[0].split('-')[0])
                amino_acid=search.group(1)
                atom=search.group(2)
                if amino_acid == 'A':
                    if atom not in {'CA','CB'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'R':
                    if atom not in {'CA','CB','CG','CD'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid in {'N','D','W','S','F','H','C'}:
                    if atom not in {'CA','CB'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid in {'Q','E'}:
                    if atom not in {'CA','CB','CG'}  or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG2','HG3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3','HG2','HG3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'G':
                    if atom not in {'CA'} or hcch_tocsy_split[0].split('-')[2] not in {'HA2','HA3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA2','HA3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'I':
                    if atom not in {'CA','CB','CG1','CG2','CD1'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB','HG12','HG13','HD1','HG2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB','HG12','HG13','HD1','HG2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'L':
                    if atom not in {'CA','CB','CG','CD1','CD2'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG','HD1','HD2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3','HG','HD1','HD2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'K':
                    if atom not in {'CA','CB','CG','CD','CE'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3','HE2','HE3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3','HG3','HG2','HD2','HD3','HE2','HE3'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'M':
                    if atom not in {'CA','CB','CG','CE'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'P':
                    if atom not in {'CA','CB','CG','CD'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB2','HB3','HG3','HG2','HE'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'T':
                    if atom not in {'CA','CB','CG2'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB','HG2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB','HG2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
                if amino_acid == 'V':
                    if atom not in {'CA','CB','CG1','CG2'} or hcch_tocsy_split[0].split('-')[2] not in {'HA','HB','HG1','HG2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} carbon or hydrogen improperly labeled')
                    if hcch_tocsy_split[0].split('-')[1] not in {'HA','HB','HG1','HG2'}:
                        print(f'Amino Acid {hcch_tocsy_split[0]} TOCSY hydrogen improperly labeled')
            except:
                print(f'Program could not analyze, stopped at \n {hcch_tocsy_lines}please check peak, correct, and rerun\n')
    print('\n')
