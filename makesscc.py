#!/usr/bin/python3
# Gruppe 03
# ProgPra WS2324
import argparse
import math
import requests
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

ATOM_TYPE_POS = (12, 16)
POSITION = (22, 26)
CHAIN_POS = 21
AA_POS = (17, 20)
X_COORD = (30, 38)
Y_COORD = (38, 46)
Z_COORD = (46, 54)
AA_ID = (22, 26)
ATOM_ID = (6, 11)
HELIX_SEQ_START = (21, 25)
HELIX_SEQ_END = (33, 37)
HELIX_CHAIN = 19
SHEET_SEQ_START = (22, 26)
SHEET_SEQ_END = (33, 37)
SHEET_CHAIN = 21

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
    ss_info = {}
    with open(id_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HELIX'):  # Get Secondary Structure-Infos for Helix
                helix_start = int(line[HELIX_SEQ_START[0]:HELIX_SEQ_START[1]].strip())
                helix_end = int(line[HELIX_SEQ_END[0]:HELIX_SEQ_END[1]].strip())
                helix_chain = line[HELIX_CHAIN]
                ss_info[(helix_chain, helix_start, helix_end)] = {'H'}
            if line.startswith('SHEET'):  # Get Secondary Structure-Infos for Sheet, else C
                sheet_start = int(line[SHEET_SEQ_START[0]:SHEET_SEQ_START[1]].strip())
                sheet_end = int(line[SHEET_SEQ_END[0]:SHEET_SEQ_END[1]].strip())
                sheet_chain = line[SHEET_CHAIN]
                ss_info[(sheet_chain, sheet_start, sheet_end)] = {'E'}
            if line.startswith('ATOM'):  # Extract Atom Info
                atom_type = line[ATOM_TYPE_POS[0]:ATOM_TYPE_POS[1]].strip()
                if atom_type == type:  # For Atoms with required Atom Type
                    chain = line[CHAIN_POS]
                    position = int(line[POSITION[0]:POSITION[1]].strip())
                    aa = line[AA_POS[0]:AA_POS[1]].strip()
                    atom_id = line[ATOM_ID[0]:ATOM_ID[1]]
                    # Get coordinates
                    x_coord = float(line[X_COORD[0]:X_COORD[1]])
                    y_coord = float(line[Y_COORD[0]:Y_COORD[1]])
                    z_coord = float(line[Z_COORD[0]:Z_COORD[1]])
                    # Create Atom_info Dict with Atom-Info-Lists as Values:
                    atom_info[atom_id] = (chain, position, atom_type, aa, x_coord, y_coord, z_coord)
            if line.startswith('MODEL'):
                counter += 1
                if counter == 2:
                    break
    max_value = atom_info[list(atom_info.keys())[-1]][1]
    return atom_info, ss_info, max_value


def calc_distance(atom_i, atom_j):
    # Read coordinates of both atoms (x1,x2,x3) (y1,y2,y3)
    i1 = atom_i[4]
    i2 = atom_i[5]
    i3 = atom_i[6]
    j1 = atom_j[4]
    j2 = atom_j[5]
    j3 = atom_j[6]
    # Return Matrixnorm of [x, y]
    return math.sqrt((j1 - i1) ** 2 + (j2 - i2) ** 2 + (j3 - i3) ** 2)


def get_sec_struct(ss_info, atom_info):
    """
    :param ss_info with keys (chain name, start_pos, end_pos)
    :param atom_info of current position
    :return: current secondary structure
    """
    for key in ss_info:
        if key[0] == atom_info[0] and key[1] <= atom_info[1] <= key[2]:
            if ss_info[key] == {'H'}:
                return 'H'
            else:
                return 'E'
    return 'C'


def calc_contacts(atom_info, distance, seq_length, ss_info):
    """
    Calculate contacts between aas based on given distance.
    :return: Dict, containing contact info for each aa
    """
    contacts = {}
    matches = {}
    distance_dict = {}
    atom_id = list(atom_info.keys())
    for i, pos_i in enumerate(atom_id):
        contacts[(pos_i,)] = {'chain': atom_info[pos_i][0],
                              'pos': atom_info[pos_i][1],
                              'serial': pos_i,
                              'aa': AA_DICT.get(atom_info[pos_i][3], 'X'),
                              'ss': get_sec_struct(ss_info, atom_info[pos_i]),
                              'global': 0,
                              'local': 0
                              }
        for j in range(0, len(atom_id)):
            if i == j:
                continue
            pos_j = atom_id[j]
            this_distance = calc_distance(atom_info[pos_i], atom_info[pos_j])
            distance_dict[atom_info[pos_i][1], atom_info[pos_j][1]] = this_distance
            if this_distance < distance:
                matches[atom_info[pos_i][1], atom_info[pos_j][1]] = 1
                if abs(atom_info[pos_i][1] - atom_info[pos_j][1]) < seq_length and atom_info[pos_i][0] == \
                        atom_info[pos_j][0]:
                    contacts[(pos_i,)]['local'] += 1
                else:
                    contacts[(pos_i,)]['global'] += 1
    return contacts, matches, distance_dict


# Download pdb file von Julius (Aufgabe get_pdb.py)
def get_pdb_file(id):
    """
    Download ptb_file from rcsb
    """
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
        header += '\n'
        file.write(header)
        for key, value in contacts.items():
            line = f"{value['chain']}\t{value['pos']}\t{value['serial']}\t{value['aa']}\t{value['ss']}\t{value['global']}\t{value['local']}"
            print(line)
            line += "\n"
            file.write(line)


def generate_contact_matrix(matches, max_value):
    matrix_size = int(max_value) + 1
    matrix = [[0] * matrix_size for _ in range(matrix_size)]

    for pos_i, pos_j in matches.keys():
        pos_i = int(pos_i)
        pos_j = int(pos_j)
        if pos_i < matrix_size and pos_j < matrix_size:
            if (pos_i, pos_j) in matches:
                matrix[pos_i][pos_j] = matches[(pos_i, pos_j)]
    return matrix


def plot_contact_matrix(contact_matrix, this_id):
    matrix = np.array(contact_matrix)
    # Calculate density of ones
    density = np.zeros_like(matrix, dtype=float)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            density[i, j] = np.sum(matrix[max(0, i - 1):min(matrix.shape[0], i + 2),
                                   max(0, j - 1):min(matrix.shape[1], j + 2)])
    # Create custom colormap from orange to dark red
    colors = plt.cm.Oranges(np.linspace(0, 1, 20))
    colors = colors[1:]  # Exclude the lightest shade of orange
    colors = np.flipud(colors)  # Reverse the colormap to start from orange
    custom_reds = mcolors.LinearSegmentedColormap.from_list('custom_reds', colors)
    # Plot the points with custom colormap
    plt.scatter(*np.where(matrix == 1)[::-1], c=density[matrix == 1], cmap=custom_reds, alpha=0.8)
    plt.gca().invert_yaxis()  # Invert y-axis to properly display the matrix

    # Add colorbar and labels
    plt.colorbar(label='Density')
    plt.xlabel('Amino-Position')
    plt.ylabel('Amino-Position')

    plt.title(f'Global and Local Contact of {this_id}s\' AminoAcids by density')

    # Show and save the plot
    plt.savefig(f"{this_id}_contact_dens_plot.png")
    print(f"{this_id}_contact_dens_plot.png successfully saved.")
    plt.show()


def generate_distance_matrix(distances):
    positions = set()
    for pos_i, pos_j in distances.keys():
        positions.add(pos_i)
        positions.add(pos_j)

    # Initialize distance matrix
    num_positions = len(positions)
    distance_matrix = np.zeros((num_positions, num_positions))

    # Fill in distance matrix
    for i, pos_i in enumerate(sorted(positions)):
        for j, pos_j in enumerate(sorted(positions)):
            if (pos_i, pos_j) in distances:
                distance_matrix[i, j] = distances[(pos_i, pos_j)]
            elif (pos_j, pos_i) in distances:
                distance_matrix[i, j] = distances[(pos_j, pos_i)]

    return distance_matrix


def plot_distance_matrix(distances, id):
    # Define custom colormap
    colors = ['white', '#FFFF33', '#FFCC33', '#FF9900', '#F3752C', '#FF0000', '#99FF00', '#990000']
    custom_cmap = mcolors.ListedColormap(colors)

    max_distance = np.max(distances)
    ticks = [0, 1, 1 + (2 * max_distance) / 9, 1 + (3 * max_distance) / 9, 1 + (4 * max_distance) / 8,
             1 + (5 * max_distance) / 8, 1 + (6 * max_distance) / 8, 1 + (7 * max_distance) / 8]
    bounds = [0, 1, 1 + (2 * max_distance) / 9, 1 + (3 * max_distance) / 9, 1 + (4 * max_distance) / 8,
              1 + (5 * max_distance) / 8, 1 + (6 * max_distance) / 8, 1 + (7 * max_distance) / 8]
    norm = mcolors.BoundaryNorm(bounds, custom_cmap.N)

    # Create heatmap of the distance matrix with custom colormap
    plt.imshow(distances, cmap=custom_cmap, norm=norm, origin='upper')

    # Add colorbar and labels
    cbar = plt.colorbar(ticks=ticks, label='Distance')
    cbar.ax.set_yticklabels(
        ['0', f'{max_distance * 1 / 8:.2f}', f'{max_distance * 2 / 8:.2f}', f'{max_distance * 3 / 8:.2f}',
         f'{max_distance * 4 / 8:.2f}', f'{max_distance * 5 / 8:.2f}', f'{max_distance * 6 / 8:.2f}',
         f'{max_distance * 7 / 8:.2f}'])
    plt.xlabel('Amino-Position')
    plt.ylabel('Amino-Position')
    plt.title(f'{id} Heatmap')

    # Show and save the plot
    plt.savefig(f"{id}_distance_plot.png")
    print(f"{id}_distance_plot.png successfully saved.")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Generate SSCC file for a given PDB ID.")
    parser.add_argument("--id", help="PDB File", required=True)
    parser.add_argument("--distance", type=float, help="Contact distance", required=True)
    parser.add_argument("--type", help="Atom type for distance calculation", required=True)
    parser.add_argument("--length", type=int, help="Sequence length for local contacts", required=True)
    parser.add_argument("--contactmatrix", help="Path to Output file for contact matrix", default=None)

    args = parser.parse_args()

    pdb_file = get_pdb_file(args.id)

    atom_info, ss_info, max_value = get_file_info(pdb_file, args.type)
    contacts, matches, distance_dict = calc_contacts(atom_info, args.distance, args.length, ss_info)
    write_sscc_table(contacts, args.id)

    if args.contactmatrix:
        contact_matrix = generate_contact_matrix(matches, max_value)
        with open(args.contactmatrix, 'w') as f:
            for line in contact_matrix:
                f.write('\t'.join(map(str, line)) + '\n')
        plot_contact_matrix(contact_matrix, args.id)

    distance_matrix = generate_distance_matrix(distance_dict)
    plot_distance_matrix(distance_matrix, args.id)


if __name__ == '__main__':
    main()
