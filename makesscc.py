#!/usr/bin/python3
import argparse
import math
import requests

ATOM_TYPE_POS = (13, 16)
POSITION = (23, 26)
CHAIN_POS = 22
AA_POS = (17, 20)
X_COORD = (31, 38)
Y_COORD = (39, 46)
Z_COORD = (47, 54)
AA_ID = (23, 26)


def get_atom_info(id_path, type):
    """
    Parse PDB File
    Extract Atom information
    :return: Dict with the positions of the atoms as keys
    """
    atom_info = {}
    with open(id_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                atom_type = line[ATOM_TYPE_POS[0]:ATOM_TYPE_POS[1]].strip()
                if atom_type == type:
                    chain = line[CHAIN_POS]
                    position = int(line[POSITION[0]:POSITION[1]].strip())
                    aa = line[AA_POS[0]:AA_POS[1]].strip()
                    # Create Atom_info Dict with Atom-Info-Lists as Values:
                    atom_info.setdefault(position, []).append((chain, position, atom_type, aa))
    return atom_info


def calc_distance(atom_i, atom_j):
    # Read coordinates of both atoms (x, y, z)
    i1 = float(atom_i[X_COORD[0]:X_COORD[1]])
    i2 = float(atom_i[Y_COORD[0]:Y_COORD[1]])
    i3 = float(atom_i[Z_COORD[0]:Z_COORD[1]])
    j1 = float(atom_j[X_COORD[0]:X_COORD[1]])
    j2 = float(atom_j[Y_COORD[0]:Y_COORD[1]])
    j3 = float(atom_j[Z_COORD[0]:Z_COORD[1]])
    # Return Matrixnorm of [i, j]
    return math.sqrt((j1 - i1) ** 2 + (j2 - i2) ** 2 + (j3 - i3) ** 2)


def calc_contacts(atom_info, distance, atom_type, seq_length):
    """
    Calculate contacts between aas based on given distance.
    :return: Dict, containing contact info for each aa
    """
    contacts = {}
    positions = list(atom_info.keys())
    for i, pos_i in enumerate(positions):
        print("Position:", pos_i)
        print("Atom Info:", atom_info)
        contacts[(pos_i,)] = {'chain': atom_info[pos_i][0],
                              'pos': pos_i,
                              'aa': atom_info[pos_i][3],
                              'ss': '',
                              'global': 0,
                              'local': 0
                              }
        for j in range(i + 1, len(positions)):
            pos_j = positions[j]
            this_distance = calc_distance(atom_info[pos_i], atom_info[pos_j])
            if this_distance <= distance:
                if abs(pos_i - pos_j) <= seq_length:
                    contacts[(pos_i,)]['local'] += 1
                else:
                    contacts[(pos_i,)]['global'] += 1
            if (pos_j,) not in contacts:
                contacts[(pos_j,)] = {'chain': atom_info[pos_j][0],
                                      'pos': pos_j,
                                      'aa': atom_info[pos_j][3],
                                      'ss': '',
                                      'global': 0,
                                      'local': 0
                                      }
            if abs(pos_i - pos_j) <= seq_length:
                contacts[(pos_j,)]['local'] += 1
            else:
                contacts[(pos_j,)]['global'] += 1
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
        for key,value in contacts.items():
            line = f"{value['chain']}\t{value['pos']}\t\t\t{value['aa']}\t{value['ss']}\t{value['global']}\t{value['local']}\n"
            print(line)
            file.write(line)


def generate_contact_matrix(contacts):
    # Find unique positions

    return 0


def main():
    parser = argparse.ArgumentParser(description="Generate SSCC file for a given PDB ID.")
    parser.add_argument("--id", help="PDB File", required=True)
    parser.add_argument("--distance", type=float, help="Contact distance", required=True)
    parser.add_argument("--type", help="Atom type for distance calculation", required=True)
    parser.add_argument("--length", type=int, help="Sequence length for local contacts", required=True)
    parser.add_argument("--contactmatrix", help="Path to Output file for contact matrix", default=None)

    args = parser.parse_args()

    pdb_file = get_pdb_file(args.id)

    atom_info = get_atom_info(pdb_file, args.type)
    contacts = calc_contacts(atom_info, args.distance, args.type, args.length)
    write_sscc_table(contacts, args.id)

    if args.contactmatrix:
        contact_matrix = generate_contact_matrix(contacts, args.contractmatrix, id)



if __name__ == '__main__':
    main()
