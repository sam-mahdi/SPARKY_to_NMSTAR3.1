import re
import sys


"""
Parameters

Sequence file should be single letter amino acid with no modifications (i.e. MSQYDGAS)

Add in peaklists, if no peaklist then leave as is. I.E. nhsqc_file = 'nhsqc.list'

save file is the name of your nmrstar hnca_file

standard_deviation_value is the cutoff for the standard deviation of yoru values, if you have a peak that is improperly labeled or has proper labeling but simply has a improper value, this will be conveyed here.
I would not recommend going above a standard deviation of 0.25
"""

sequence_file='seq.txt'

nhsqc_file='G_ML_Nhsqc.list'
hnca_file='G_ML_HNCA.list'
hncacb_file='G_ML_HNCACB.list'
hncoca_file='G_ML_HNCAi-1.list'
hnco_file='G_ML_HNCO.list'
hbhaconh_file='G_ML_HBHACONH.list'
chsqc_file='G_ML_Chsqc.list'
cch_tocsy_file='G_ML_CCH_TOCSY.list'
hcch_tocsy_file='G_ML_HCCH_TOCSY.list'

save_file='G_ML_STAR3.1.txt'

standard_deviation_value=0.25

print('Starting Program')

def sequence_list():
    seq_list=[]
    counter=0
    with open(sequence_file) as file:
        for values in file:
            strip_values=values.strip()
            for entries in strip_values:
                counter+=1
                seq_list.append(f'{counter}{entries}')
    return seq_list
def regex_list():
    new_list=[]
    Alanine_Values=['C','CA','CB','HA','HB','HN','N']
    Arganine_Values=['C','CA','CB','CD','CG','HA','HB2','HB3','HD2','HD3','HG2','HG3','HN','N']
    Aspartic_Acid_Values=['C','CA','CB','HA','HB2','HB3','HN','N']
    Glutamine_Values=['C','CA','CB','CG','HA','HB2','HB3','HE21','HE22','HG2','HG3','HN','N','NE2']
    Glycine_Values=['C','CA','HA2','HA3','HN','N']
    Isoleucine_Values=['C','CA','CB','CD1','CG1','CG2','HA','HB','HD1','HG12','HG13','HG2','HN','N']
    Luecine_Values=['C','CA','CB','CD1','CD2','CG','HA','HB2','HB3','HD1','HD2','HG','HN','N']
    Lysine_Values=['C','CA','CB','CD','CE','CG','HA','HB2','HB3','HD2','HD3','HE2','HE3','HG2','HG3','HN','N']
    Methionine_Values=['C','CA','CB','CE','CG','HA','HB2','HB3','HE','HG2','HG3','HN','N']
    Proline_Values=['C','CA','CB','CD','CG','HA','HB2','HB3','HD2','HD3','HG2','HG3','HN','N']
    Serine_Values=['C','CA','CB','HA','HB2','HB3','HN','N']
    Threanine_Values=['C','CA','CB','CG2','HA','HB','HG2','HN','N']
    Valine_Values=['C','CA','CB','CG1','CG2','HA','HB','HG2','HN','N']
    for amino_acids in sequence_list():
        if amino_acids[-1] == 'M':
            methionine=[amino_acids+'-'+atom for atom in Methionine_Values]
            new_list+=methionine
        if amino_acids[-1] in {'D','N','C','F','S','W','Y','H'}:
            aspartic_acid=[amino_acids+'-'+atom for atom in Aspartic_Acid_Values]
            new_list+=aspartic_acid
        if amino_acids[-1] == 'A':
            alanine=[amino_acids+'-'+atom for atom in Alanine_Values]
            new_list+=alanine
        if amino_acids[-1] in {'Q','E'}:
            glutamine=[amino_acids+'-'+atom for atom in Glutamine_Values]
            new_list+=glutamine
        if amino_acids[-1] == 'G':
            glycine=[amino_acids+'-'+atom for atom in Glycine_Values]
            new_list+=glycine
        if amino_acids[-1] == 'R':
            arginine=[amino_acids+'-'+atom for atom in Arganine_Values]
            new_list+=arginine
        if amino_acids[-1] == 'I':
            isoleucine=[amino_acids+'-'+atom for atom in Isoleucine_Values]
            new_list+=isoleucine
        if amino_acids[-1] == 'L':
            leucine=[amino_acids+'-'+atom for atom in Luecine_Values]
            new_list+=leucine
        if amino_acids[-1] == 'K':
            lysine=[amino_acids+'-'+atom for atom in Lysine_Values]
            new_list+=lysine
        if amino_acids[-1] == 'P':
            proline=[amino_acids+'-'+atom for atom in Proline_Values]
            new_list+=proline
        if amino_acids[-1] == 'T':
            threanine=[amino_acids+'-'+atom for atom in Threanine_Values]
            new_list+=threanine
        if amino_acids[-1] == 'V':
            valine=[amino_acids+'-'+atom for atom in Valine_Values]
            new_list+=valine
    return new_list

def NHSQC_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(nhsqc_file) as nhsqc:
        for nhsqc_lines in nhsqc:
            if re.search('^\w\d+\w+',nhsqc_lines.strip()) is None:
                continue
            nhsqc_split=nhsqc_lines.strip().split()
            if amino_acid+residue_number+atom == nhsqc_split[0].split('-')[0]:
                temp_list.append(float(nhsqc_split[1]))
                spectra_list.append('NHSQC')
            if atom == nhsqc_split[0].split('-')[1] and amino_acid+residue_number == (re.search('[A-Z]\d+',nhsqc_split[0].split('-')[0])).group(0):
                temp_list.append(float(nhsqc_split[2]))
                spectra_list.append('NHSQC')
def HNCA_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(hnca_file) as hnca:
        for hnca_lines in hnca:
            if re.search('^\w\d+\w+',hnca_lines.strip()) is None:
                continue
            hnca_split=hnca_lines.strip().split()
            if amino_acid+residue_number == hnca_split[0].split('-')[0][0:-1] and atom == hnca_split[0].split('-')[1]:
                temp_list.append(float(hnca_split[2]))
                spectra_list.append('HNCA')
            if amino_acid+residue_number+atom == hnca_split[0].split('-')[1]:
                temp_list.append(float(hnca_split[2]))
                spectra_list.append('HNCA')
def HNCACB_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(hncacb_file) as hncacb:
        for hncacb_lines in hncacb:
            if re.search('^\w\d+\w+',hncacb_lines.strip()) is None:
                continue
            hncacb_split=hncacb_lines.strip().split()
            if amino_acid+residue_number == hncacb_split[0].split('-')[0][0:-1] and atom == hncacb_split[0].split('-')[1]:
                temp_list.append(float(hncacb_split[2]))
                spectra_list.append('HNCACB')
            if amino_acid+residue_number+atom == hncacb_split[0].split('-')[1]:
                temp_list.append(float(hncacb_split[2]))
                spectra_list.append('HNCACB')
def HNCOCA_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(hncoca_file) as hnca1:
        for hnca1_lines in hnca1:
            if re.search('^\w\d+\w+',hnca1_lines.strip()) is None:
                continue
            hnca1_split=hnca1_lines.strip().split()
            if amino_acid+residue_number+atom == hnca1_split[0].split('-')[1]:
                temp_list.append(float(hnca1_split[2]))
                spectra_list.append('HNCOCA')
def HNCO_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(hnco_file) as hnco:
        for hnco_lines in hnco:
            if re.search('^\w\d+\w+',hnco_lines.strip()) is None:
                continue
            hnco_split=hnco_lines.strip().split()
            if amino_acid+residue_number+atom == hnco_split[0].split('-')[1]:
                temp_list.append(float(hnco_split[2]))
                spectra_list.append('HNCO')
def CHSQC_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(chsqc_file) as chsqc:
        for chsqc_lines in chsqc:
            if re.search('^\w\d+\w+',chsqc_lines.strip()) is None:
                continue
            chsqc_split=chsqc_lines.strip().split()
            if amino_acid+residue_number+atom == chsqc_split[0].split('-')[0]:
                temp_list.append(float(chsqc_split[1]))
                spectra_list.append('CHSQC')
            if atom == chsqc_split[0].split('-')[1] and amino_acid+residue_number == (re.search('[A-Z]\d+',chsqc_split[0].split('-')[0])).group(0):
                temp_list.append(float(chsqc_split[2]))
                spectra_list.append('CHSQC')
def HBHACONH_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(hbhaconh_file) as hbhaconh:
        for hbhaconh_lines in hbhaconh:
            if re.search('^\w\d+\w+',hbhaconh_lines.strip()) is None:
                continue
            hbhaconh_split=hbhaconh_lines.strip().split()
            if amino_acid+residue_number+atom == hbhaconh_split[0].split('-')[1]:
                temp_list.append(float(hbhaconh_split[2]))
                spectra_list.append('HBHACONH')
def CCH_TOCSY_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(cch_tocsy_file) as ccc_tocsy:
        for ccc_tocsy_lines in ccc_tocsy:
            if re.search('^\w\d+\w+',ccc_tocsy_lines.strip()) is None:
                continue
            ccc_tocsy_split=ccc_tocsy_lines.strip().split()
            if amino_acid+residue_number == (re.search('[A-Z]\d+',ccc_tocsy_split[0].split('-')[0])).group(0) and atom == ccc_tocsy_split[0].split('-')[1]:
                temp_list.append(float(ccc_tocsy_split[2]))
                spectra_list.append('CCH TOCSY')
def HCCH_TOCSY_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list):
    with open(hcch_tocsy_file) as hcch_tocsy:
        for hcch_tocsy_lines in hcch_tocsy:
            if re.search('^\w\d+\w+',hcch_tocsy_lines.strip()) is None:
                continue
            hcch_tocsy_split=hcch_tocsy_lines.strip().split()
            if amino_acid+residue_number == (re.search('[A-Z]\d+',hcch_tocsy_split[0].split('-')[0])).group(0) and atom == hcch_tocsy_split[0].split('-')[2]:
                temp_list.append(float(hcch_tocsy_split[3]))
                spectra_list.append('HCCH TOCSY')

def check_peaklist_labels():
    print('Checking proper labeling and formatting in peaklists')
    import checker as ch
    if nhsqc_file != '':
        ch.NHSQC_checker(nhsqc_file)
    if hnca_file != '':
        ch.HNCA_checker(hnca_file)
    if hncacb_file != '':
        ch.HNCACB_checker(hncacb_file)
    if hncoca_file != '':
        ch.HNCOCA_checker(hncoca_file)
    if hnco_file != '':
        ch.HNCO_checker(hnco_file)
    if chsqc_file != '':
        ch.CHSQC_checker(chsqc_file)
    if hbhaconh_file != '':
        ch.HBHACONH_checker(hbhaconh_file)
    if cch_tocsy_file != '':
        ch.CCH_TOCSY_checker(cch_tocsy_file)
    if hcch_tocsy_file != '':
        ch.HCCH_TOCSY_checker(hcch_tocsy_file)


def nmrstar_file():
    print('Generating NMRSTAR File (Takes a minute depending on size of protein)')
    aa_dict={'D':'Asp','T':'Thr','S':'Ser','E':'Glu','P':'Pro','G':'Gly','A':'Ala','C':'Cys','V':'Val',
    'M':'Met','I':'Ile','L':'Leu','Y':'Tyr','F':'Phe','H':'His','K':'Lys','R':'Arg','W':'Trp','Q':'Gln','N':'Asn'}
    temp_list=[]
    nmrstar_list=[]
    spectra_list=[]
    counter=0
    for values in regex_list():
        split_values=values.split('-')
        amino_acid=split_values[0][-1]
        residue_number=split_values[0][0:-1]
        atom=split_values[1]
        if nhsqc_file != '':
            NHSQC_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if hnca_file != '':
            HNCA_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if hncacb_file != '':
            HNCACB_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if hncoca_file != '':
            HNCOCA_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if hnco_file != '':
            HNCO_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if chsqc_file != '':
            CHSQC_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if hbhaconh_file != '':
            HBHACONH_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if cch_tocsy_file != '':
            CCH_TOCSY_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)
        if hcch_tocsy_file != '':
            HCCH_TOCSY_peaklist(amino_acid,residue_number,atom,temp_list,spectra_list)

        if len(temp_list) == 0:
            continue
        counter+=1
        average=(sum(temp_list))/(len(temp_list))
        standard_deviation = (sum([((x - average) ** 2) for x in temp_list]) / len(temp_list))**0.5
        if standard_deviation > standard_deviation_value:
            print(values)
            print(spectra_list)
            print(temp_list)
        if atom[0] == 'C':
            nmrstar_list.append(f'{counter} {residue_number} {residue_number} {aa_dict[amino_acid]} {atom} C 13 {round(average,3)} {round(standard_deviation,2)} 0 1')
        if atom[0] == 'N':
            nmrstar_list.append(f'{counter} {residue_number} {residue_number} {aa_dict[amino_acid]} {atom} N 15 {round(average,3)} {round(standard_deviation,2)} 0 1')
        if atom[0] == 'H':
            nmrstar_list.append(f'{counter} {residue_number} {residue_number} {aa_dict[amino_acid]} {atom} H 1 {round(average,3)} {round(standard_deviation,2)} 0 1')
        temp_list.clear()
        spectra_list.clear()
    rows = [line.split() for line in nmrstar_list]
    columns = zip(*rows)
    col_widths = [max(len(e) for e in column) for column in columns]

    new_data = []
    for row in rows:
        new_row = ''
        for e, width in zip(row, col_widths):
            new_row += f"{e:>{width}} "
        new_data.append(new_row)

    with open(save_file,'w') as file:
        file.write('loop_\n    _Atom_chem_shift.ID\n    _Atom_chem_shift.Comp_index_ID\n    _Atom_chem_shift.Seq_ID\n    _Atom_chem_shift.Comp_ID\n    _Atom_chem_shift.Atom_ID\n    _Atom_chem_shift.Atom_type\n    _Atom_chem_shift.Atom_isotope_number\n    _Atom_chem_shift.Val\n    _Atom_chem_shift.Val_err\n    _Atom_chem_shift.Ambiguity_code\n     _Atom_chem_shift.Assigned_chem_shift_list_ID\n')
        for row in new_data:
            file.write(f'{row}\n')
        file.write('stop_')




def compare_to_bmrb():
    print('Comparing to BMRB Values')
    with open(save_file) as file:
        for lines in file:
            if re.search('^\d+',lines.strip()) is None:
                continue
            star_residue=lines.split()[1]
            star_amino_acid=lines.split()[3].upper()
            star_atom=lines.split()[4]
            star_value=lines.split()[7]
            with open('bmrb.csv') as file2:
                for lines2 in file2:
                    if lines2 == '\n':
                        continue
                    split_lines=lines2.split(',')
                    residue=split_lines[0]
                    if residue == 'comp_id':
                        continue
                    atom=split_lines[1]
                    chemical_shift=float(split_lines[5])
                    std=float(split_lines[6])
                    lower_half=chemical_shift-std
                    upper_half=chemical_shift+std
                    if star_amino_acid == residue and star_atom == atom:
                        if float(star_value) > (upper_half+(0.05*upper_half)):
                            print(f'{star_amino_acid} {star_residue} {star_atom} value {star_value} is outside of range {round(lower_half,3)}-{round(upper_half,3)}\n')
                        if float(star_value) < (lower_half-(0.05*lower_half)):
                            print(f'{star_amino_acid} {star_residue} {star_atom} value {star_value} is outside of range {round(lower_half,3)}-{round(upper_half,3)}\n')

def main_loop():
    check_peaklist_labels()
    print('Check Complete')
    pause=input('If no errors in formatting, please type continue. \nOtherwise, type quit, correct the errors, and rerun: ')
    if pause == 'quit':
        sys.exit()
    nmrstar_file()
    compare_to_bmrb()
main_loop()
