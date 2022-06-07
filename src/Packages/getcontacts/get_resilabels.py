#!/usr/bin/env python
"""
Take .txt output of multiple structure alignment between PDBs as generated by gesamt
(http://www.ccp4.ac.uk/dist/html/gesamt.html) and output label files to use with GetContacts.
"""

import os


class Aligned_Residues:
    def __init__(self, processed_line, include_nonaligned=False):
        resnames = [resname for protein, ss, resname in processed_line]
        self.ss_colormap = {"H": "red", "S": "yellow"}
        if include_nonaligned:
            self.to_print = True
        else:
            self.to_print = sum([bool(resname) for resname in resnames]) > 1
        print(self.to_print)
        self.generic = "-".join(resnames)
        self.protein_map = {
            protein: {"resname": resname, "ss": ss, "to_print": (len(resname) > 0)}
            for protein, ss, resname in processed_line
        }

    def get_line(self, protein):
        if not self.to_print or not self.protein_map[protein]["to_print"]:
            return ""
        line = "{}\t{}".format(self.protein_map[protein]["resname"], self.generic)
        if (
            self.protein_map[protein]["ss"]
            and self.protein_map[protein]["ss"] in self.ss_colormap
        ):
            line += "\t{}".format(self.ss_colormap[self.protein_map[protein]["ss"]])
        line += "\n"
        return line


def parse_two_queries(gesamt_lines, include_nonaligned, proteins=None):
    if proteins is None:
        proteins = []
        for idx, line in enumerate(gesamt_lines):
            if (
                "reading QUERY structure" in line
                or "reading TARGET structure" in line
                or "reading FIXED structure" in line
                or "reading MOVING structure" in line
            ):
                protein_filename = line.split("'")[1]
                protein, _ = os.path.splitext(protein_filename)
                proteins.append(protein)

    # get the index of the first line of the alignment and the index of the last line of the alignment
    alignment_start_idx = None
    alignment_end_idx = None
    for idx, line in enumerate(gesamt_lines):
        if not alignment_start_idx:
            line_tokens = [token.strip() for token in line.split("|")]
            # the header line looks like this: "|    Query    |  Dist.(A)  |   Target    |"
            if (
                ("Query" in line_tokens or "FIXED" in line_tokens)
                and "Dist.(A)" in line_tokens
                and ("Target" in line_tokens or "MOVING" in line_tokens)
            ):
                alignment_start_idx = idx + 2
        else:
            if "'" in line:
                alignment_end_idx = idx
                break

    # grab all lines in the file belonging to the alignment
    alignment_lines = gesamt_lines[alignment_start_idx:alignment_end_idx]

    # slice out three tokens for each line, corresponding to Query Residue, Distance, and Target Residue
    # for example: "|H- A:LEU  75 | <**0.82**> |H- A:LEU  65 |" -> ["H- A:LEU  75", "H- A:LEU  65"]
    alignment_lines = [
        [line.split("|")[idx] for idx in [1, 3]] for line in alignment_lines
    ]

    aligned_residues = []
    for line in alignment_lines:
        new_line = []
        for protein, token in zip(proteins, line):
            resname = (
                token[-10:-5] + ":" + str(int(token[-5:])) if token.strip() else ""
            )
            ss = token[0].strip()
            new_line.append((protein, ss, resname))
        print(new_line)
        aligned_residues.append(Aligned_Residues(new_line, include_nonaligned))

    return aligned_residues


def parse_more_than_two_queries(gesamt_lines, include_nonaligned, proteins=None):
    if proteins is None:
        proteins = []
        for idx, line in enumerate(gesamt_lines):
            if "... reading file" in line:
                protein_filename = line.split("'")[1]
                protein, _ = os.path.splitext(protein_filename)
                proteins.append(protein)

    alignment_start_idx = None
    alignment_end_idx = None
    for idx, line in enumerate(gesamt_lines):
        if not alignment_start_idx:
            line_tokens = [token.strip() for token in line.split("|")]
            if "Disp." in line_tokens:
                alignment_start_idx = idx + 2
        else:
            if "'" in line:
                alignment_end_idx = idx
                break

    # grab all lines in the file belonging to the alignment
    alignment_lines = gesamt_lines[alignment_start_idx:alignment_end_idx]

    aligned_residues = []
    for line in alignment_lines:
        # first, replace any instance of the weird delimiter with the asterisk in the center: "6.034 |*|  A:CYS 341 |*|  A:MET 456 |*|  D:LEU 559" -> "6.034 | |  A:CYS 341 | |  A:MET 456 | |  D:LEU 559"
        line = line.replace("|*|", "| |")

        # next, remove the line delimeters: "6.034 | |  A:CYS 341 | |  A:MET 456 | |  D:LEU 559" -> "['6.034', 'A:CYS 341', 'A:MET 456', 'D:LEU 559']"
        line = [token.strip() for token in line.split("| |")]

        # remove the first (displacement) column from each line: "['6.034', 'A:CYS 341', 'A:MET 456', 'D:LEU 559']" -> "['A:CYS 341', 'A:MET 456', 'D:LEU 559']"
        line = line[1:]

        new_line = []
        for protein, token in zip(proteins, line):
            if len(token.split("|")) == 1:
                ss = ""
                resname = token
            elif len(token.split("|")) == 2:
                ss, resname = token.split("|")
            resname = resname[:5] + ":" + str(int(resname[5:])) if resname else ""
            new_line.append((protein, ss, resname))
        print(new_line)
        aligned_residues.append(Aligned_Residues(new_line))

    return aligned_residues


def main(argv=None):
    # Set up command line arguments
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input_gesamt",
        required=True,
        type=str,
        help="A .txt file produced by redirecting output from gesamt to a file",
    )
    parser.add_argument(
        "--output_path", required=True, type=str, help="Path to output directory"
    )
    parser.add_argument(
        "--proteins",
        required=False,
        type=str,
        nargs="+",
        help='Optionally, provide a space-delimited list of names of all proteins in the gesamt input file. Protein names should be specified in the order that the proteins were inputted to the gesamt program. For example, if you called gesamt with "gesamt foo1.pdb foo2.pdb foo3.pdb", you can specify "--proteins foo1 foo2 foo3"',
    )
    parser.add_argument(
        "--include_nonaligned",
        required=False,
        dest="include_nonaligned",
        action="store_true",
        help="Optionally, include in the label file residues that are not aligned to any other residue. This option is not recommended.",
    )
    parser.set_defaults(include_nonaligned=False)

    args = parser.parse_args(argv)
    input_gesamt = args.input_gesamt
    output_path = args.output_path
    proteins = args.proteins
    include_nonaligned = args.include_nonaligned

    with open(input_gesamt, "r") as r_open:
        gesamt_lines = [
            line.strip() for line in r_open.readlines() if len(line.strip()) > 0
        ]

    # parse the input file differently depending on whether there are 2 queries or >2 queries
    gesamt_text = "".join(gesamt_lines)
    if (
        "reading QUERY structure" in gesamt_text
        or "reading FIXED structure" in gesamt_text
    ):
        aligned_residues = parse_two_queries(
            gesamt_lines, include_nonaligned, proteins=proteins
        )
    else:
        aligned_residues = parse_more_than_two_queries(
            gesamt_lines, include_nonaligned, proteins=proteins
        )

    os.system("mkdir -p {}".format(output_path))

    for protein in proteins:
        # write output file
        with open("{}/{}.label".format(output_path, protein), "w+") as w_open:
            for residue_set in aligned_residues:
                w_open.write(residue_set.get_line(protein))


if __name__ == "__main__":
    main()