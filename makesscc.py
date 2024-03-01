#!/usr/bin/python3
import argparse
import math
import requests

# AA 3 -> 1, glob, local, helix

ATOM_TYPE_POS = (13, 16)
POSITION = (23, 26)
CHAIN_POS = 22
AA_POS = (18, 20)
X_COORD = (31, 38)
Y_COORD = (39, 46)
Z_COORD = (47, 54)
AA_ID = (23, 26)
ATOM_ID = (7, 11)
HELIX_SEQ_START = (22, 25)
HELIX_SEQ_END = (34, 37)
HELIX_CHAIN = 20
SHEET_SEQ_START = (23, 26)
SHEET_SEQ_END = (34, 37)
SHEET_CHAIN = 22

AA_DICT = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


def get_file_info(id_path, type):
    """
    Parse PDB File
    Extract Atom information
    :return: Dict with the positions of the atoms as keys
    """
    counter = 0
    atom_info = {}
    ss_info = {'helix': [], 'sheet': []}
    with open(id_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HELIX'):  # Get Secondary Structure-Infos for Helix
                helix_start = int(line[HELIX_SEQ_START[0] - 1:HELIX_SEQ_START[1]].strip())
                helix_end = int(line[HELIX_SEQ_END[0] - 1:HELIX_SEQ_END[1]].strip())
                belonging_chain = line[HELIX_CHAIN - 1]
                ss_info['helix'].append([helix_start, helix_end, belonging_chain])
            if line.startswith('SHEET'):  # Get Secondary Structure-Infos for Sheet, else C
                sheet_start = int(line[SHEET_SEQ_START[0] - 1:SHEET_SEQ_START[1]].strip())
                sheet_end = int(line[SHEET_SEQ_END[0] - 1:SHEET_SEQ_END[1]].strip())
                belonging_chain = line[SHEET_CHAIN - 1]
                ss_info['sheet'].append([sheet_start, sheet_end, belonging_chain])
            if line.startswith('ATOM'):  # Extract Atom Info
                atom_type = line[ATOM_TYPE_POS[0] - 1:ATOM_TYPE_POS[1]].strip()
                if atom_type == type:  # For Atoms with required Atom Type
                    chain = line[CHAIN_POS - 1]
                    position = int(line[POSITION[0] - 1:POSITION[1]].strip())
                    aa = line[AA_POS[0] - 1:AA_POS[1]].strip()
                    atom_id = line[ATOM_ID[0] - 1:ATOM_ID[1]]
                    # Get coordinates
                    x_coord = float(line[X_COORD[0] - 1:X_COORD[1]])
                    y_coord = float(line[Y_COORD[0] - 1:Y_COORD[1]])
                    z_coord = float(line[Z_COORD[0] - 1:Z_COORD[1]])
                    print(x_coord, y_coord)
                    # Create Atom_info Dict with Atom-Info-Lists as Values:
                    atom_info[atom_id] = (chain, position, atom_type, aa, x_coord, y_coord, z_coord)
            if line.startswith('MODEL'):
                counter += 1
                if counter == 2:
                    break
    print(ss_info)
    return atom_info, ss_info


def calc_distance(atom_i, atom_j):
    # Read coordinates of both atoms (x, y, z)
    i1 = atom_i[4]
    i2 = atom_i[5]
    i3 = atom_i[6]
    j1 = atom_j[4]
    j2 = atom_j[5]
    j3 = atom_j[6]
    # Return Matrixnorm of [i, j]
    return math.sqrt((j1 - i1) ** 2 + (j2 - i2) ** 2 + (j3 - i3) ** 2)


def get_secondary_structure(pos, chain, ss_info):
    print(ss_info)
    helix_info = ss_info['helix']
    sheet_info = ss_info['sheet']

    for info in helix_info:
        if info[0] <= int(pos) <= info[1] and info[2] == chain:
            return 'H'
    for info in sheet_info:
        if info[0] <= int(pos) <= info[1] and info[2] == chain:
            return 'E'
    return 'C'  # Coil


def calc_contacts(atom_info, distance, seq_length, ss_info):
    """
    Calculate contacts between aas based on given distance.
    :return: Dict, containing contact info for each aa
    """
    contacts = {}
    atom_id = list(atom_info.keys())
    for i, pos_i in enumerate(atom_id):
        contacts[(pos_i,)] = {'chain': atom_info[pos_i][0],
                              'pos': atom_info[pos_i][1],
                              'serial': pos_i,
                              'aa': AA_DICT[atom_info[pos_i][3]],
                              'ss': get_secondary_structure(pos_i, atom_info[pos_i][0], ss_info),
                              'global': 0,
                              'local': 0
                              }
        for j in range(0, len(atom_id)):
            if i == j:
                continue
            pos_j = atom_id[j]
            this_distance = calc_distance(atom_info[pos_i], atom_info[pos_j])
            if this_distance < distance:
                if abs(atom_info[pos_i][1] - atom_info[pos_j][1]) < seq_length and atom_info[pos_i][0] == atom_info[pos_j][0]:
                    contacts[(pos_i,)]['local'] += 1
                else:
                    contacts[(pos_i,)]['global'] += 1
    return contacts


# Download pdb file von Julius (Aufgabe get_pdb.py)
def get_pdb_file(id):
    url = f'https://files.rcsb.org/download/{id}.pdb'
    response = requests.get(url)
    file_path = f"{id}.pdb"
    with open(file_path, 'w') as file:
        file.write(response.text)
    return file_path


def write_sscc_table(contacts, id):
    output_file = f'{id}.sscc'
    with open(output_file, 'w') as file:
        header = "chain\tpos\tserial\taa\tss\tglobal\tlocal"
        print(header)
        file.write(header)
        for key, value in contacts.items():
            line = f"{value['chain']}\t{value['pos']}\t{value['serial']}\t{value['aa']}\t{value['ss']}\t{value['global']}\t{value['local']}"
            print(line)
            file.write(line)


def generate_contact_matrix():

    return


def main():
    parser = argparse.ArgumentParser(description="Generate SSCC file for a given PDB ID.")
    parser.add_argument("--id", help="PDB File", required=True)
    parser.add_argument("--distance", type=float, help="Contact distance", required=True)
    parser.add_argument("--type", help="Atom type for distance calculation", required=True)
    parser.add_argument("--length", type=int, help="Sequence length for local contacts", required=True)
    parser.add_argument("--contactmatrix", help="Path to Output file for contact matrix", default=None)

    args = parser.parse_args()

    pdb_file = get_pdb_file(args.id)

    atom_info, ss_info = get_file_info(pdb_file, args.type)
    contacts = calc_contacts(atom_info, args.distance, args.length, ss_info)
    write_sscc_table(contacts, args.id)

    if args.contactmatrix:
        contact_matrix = generate_contact_matrix()
        with open(args.contactmatrix, 'w') as f:
            for line in contact_matrix:
                f.write('\t'.join(map(str, line)) + '\n')


if __name__ == '__main__':
    main()
