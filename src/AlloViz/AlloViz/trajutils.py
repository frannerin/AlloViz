"""Module containing functions for trajectory processing

Functions in this module are used by :class:`AlloViz.Protein` for input file processing
or by other functions herein to process and prepare the passed input files for analysis
with the rest of AlloViz's functionalities.

"""

import os

import requests, parmed
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisFromFunction
import numpy as np
from multiprocess import Pool

from . import utils


def standardize_resnames(protein, **kwargs):
    r"""Change the residue names of a protein to standard 3-letter codes

    Change the residue names of the passed :class:`~MDAnalysis.core.universe.Universe` to
    standard 3-letter names and also put them in the same segment (the first segment of
    the Universe by default) to avoid problems with packages not recognizing non-standard
    residue names or producing unexpected results if residues are in various segments.

    :mod:`Bio.SeqUtils` is used to change the codes and also retrieve the standard 3-to-1
    and 1-to-3 residue code mapping, which is also extended with
    :data:`MDAnalysis.lib.util.inverse_aa_codes`.

    Parameters
    ----------
    protein : :class:`~MDAnalysis.core.groups.AtomGroup`
        Atoms of the protein for which to standardize residue names.
    **kwargs
        `special_res` can be passed as an optional kwarg, and it should be a dictionary
        containing a mapping of special residue 3/4-letter code(s) present in the
        structure to the corresponding standard 1-letter code(s).
    """
    from Bio.SeqUtils import seq1, seq3
    from Bio.SeqUtils import IUPACData
    from MDAnalysis.lib.util import inverse_aa_codes

    # Make a dictionary of the 3-to-1 letter residue codes mapping
    # It is made with a combination of the info from Bio.SeqUtils.IUPACData and MDAnalysis.lib.util.inverse_aa_codes and, if applicable, the special_res kwarg dict
    res_d = {}
    res_d.update(IUPACData.protein_letters_3to1_extended)
    res_d.update(inverse_aa_codes)
    if "special_res" in kwargs and isinstance(kwargs["special_res"], dict):
        res_d.update(kwargs["special_res"])

    # Set all residues' segment to the same one to avoid PyInteraph2_Energy problems
    protein.residues.segments = protein.segments[0]

    # Try to standardize all residues' names and, if one fails, raise an Exception to ask for the 3-to-1 mapping
    try:
        for res in protein.residues:
            res.resname = seq3(seq1(res.resname, custom_map=res_d)).upper()
    except:
        raise Exception(
            f"""
resname(s) {[res for res in set(protein.residues.resnames) if res not in res_d]} could not be translated to a standard name.
Please provide a mapping from 3/4-letter code to 1-letter code as a dictionary with keyword argument 'special_res' upon initialization.
              """
        )


def get_GPCRdb_numbering(protein):
    r"""Retrieve the GPCRdb-assigned generic numbering of a GPCR structure

    A PDB file with GPCR generic numbering of residues in the B-factor column of the file
    is retrieved through the `GPCRdb <https://gpcrdb.org/>`_ API and the GPCRdb-scheme 
    generic numbers are transfered to the `protein` CA atoms' `tempfactors`.

    GPCRdb assigns Ballesteros-Weinstein generic numbers to the N atoms' B-factors, and
    the GPCRdb-scheme ones to the CA atoms'. GPCRdb "bulges" are marked by repeating the
    generic number but making it negative (negative B-factor instead of a three-decimal
    B-factor) and must be transformed by `-b + 0.001`. Some residual `1.00` B-factors
    remain and must be taken care of.

    Parameters
    ----------
    protein : :class:`~MDAnalysis.core.groups.AtomGroup`
        Atoms of the protein for which to retrieve generic numbers and that will be
        transformed (its CA atoms' `tempfactors`) in-place.

    Notes
    -----
    GPCRdb provides a PyMOL file to expedite the visualization of the two schemes of
    generic numbers as text tags over the structure, which is where the information about
    what is stored in which atom's B-factor column is taken from:

    .. code-block::
    
        print "GPCRdb script labeling generic numbers on the annotated pdb structure\nKey bindings for labels:\nF1 - show generic numbers\nF2 - show Ballesteros-Weinstein numbers\nF3 - clear labels"
        label n. CA or n. N
        cmd.set_key('F1', 'label n. CA & (b >-8.1 and  b <8.1), str("%1.2f" %b).replace(".","x") if b > 0 else str("%1.3f" %(-b + 0.001)).replace(".", "x")')
        cmd.set_key('F2', 'label n. N & (b > 0 and  b <8.1), "%1.2f" %b')
        cmd.set_key('F3', 'label n. CA or n. N')
    """
    from io import StringIO as IO

    ### USING CONTEXT MANAGERS
    # This has to be a function itself instead of linear code because there are problems with MDAnalysis' NamedStreams not closing properly or something and it throws a Warning
    def get_GPCRdb_result(protein):
        # An MDAnalysis' NamedStream is used to write the passed `protein` as a PDB file and pass it to the GPCRdb API; a PDB with the generic numbers is retrieved directly as a response
        with IO("") as protf, mda.lib.util.NamedStream(protf, "prot.pdb") as f:
            protein.write(f, file_format="pdb")
            GPCRdb_output = requests.post(
                "https://gpcrdb.org/services/structure/assign_generic_numbers",
                files={"pdb_file": protf.getvalue()},
            )
            return GPCRdb_output.text

    # Another MDAnalysis' NamedStream is used to be able to use the obtained generic-numbered PDB as a file from which a Universe can be constructed to extract the B-factors
    out = get_GPCRdb_result(protein)
    with mda.lib.util.NamedStream(IO(out), "gen_numb.pdb") as gen_numb:
        # GPCRdb-scheme generic numbers are assigned to the CA atoms' and residual `1.00` B-factors must be taken care of
        # GPCRdb "bulges" are marked by repeating the generic number but making it negative (negative B-factor instead of a three-decimal B-factor) and must be transformed by `-b + 0.001`
        nums = mda.Universe(gen_numb).select_atoms("name CA").tempfactors
        GPCRdb_nums = np.array(
            [
                round(bfactor, 2)
                if bfactor > 1
                else round(-bfactor + 0.001, 3)
                if bfactor < 0
                else 0
                for bfactor in nums
            ]
        )
        # Transform the passed `protein` in-place
        protein.select_atoms("name CA").tempfactors = GPCRdb_nums

    return

    ### WITHOUT CONTEXT MANAGERS


#     protf = IO()
#     f = mda.lib.util.NamedStream(protf, "prot.pdb")

#     protein.write(f, file_format="pdb")
#     GPCRdb_output = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers',
#                                   files = {'pdb_file': protf.getvalue()})
#     out = GPCRdb_output.text

#     gen_numb = mda.lib.util.NamedStream(IO(out), "gen_numb.pdb")
#     nums = mda.Universe(gen_numb).select_atoms("name CA").tempfactors
#     GPCRdb_nums = np.array([round(bfactor, 2) if bfactor>1 else round(-bfactor+0.001, 3) if bfactor<0 else 0 for bfactor in nums])
#     protein.select_atoms("name CA").tempfactors = GPCRdb_nums


def download_GPCRmd_files(GPCRmdid, path):
    r"""Download the files of a GPCRmd dynid stored in the database

    The website of the corresponding `GPCRmd <https://gpcrmd.org/>`_ ID is scanned to
    retrieve the downloadable files and they are downloaded in parallel in the passed
    path with the same name that they have in the database.

    Parameters
    ----------
    GPCRmdid : int
        Dynid from the GPCRmd database.
    path : str
        Path, relative or absolute, in which to download the files.

    Notes
    -----
    The force-field parameters file is downloaded as a tar file and is extracted and also
    transformed in-place to avoid downstream problems with the gRINN construction method.
    """
    import sys, re
    import tarfile, fileinput
    from bs4 import BeautifulSoup
    import urllib.request as pwget

    web = "https://submission.gpcrmd.org"

    # Establish the URL from which the files can be accessed and find the donwload links for the file extensions of interest
    html = requests.get(f"{web}/dynadb/dynamics/id/{GPCRmdid}/")
    soup = BeautifulSoup(html.text, features="html.parser").find_all("a")
    links = [
        link.get("href")
        for link in soup
        if re.search("(xtc|pdb|psf|prm)", link.get("href"))
    ]

    # Download the files in parallel and save them in `path` with the same original name that they have in GPCRmd (last part of `link`)
    mypool = Pool()

    get_name = lambda link: f"{path}/{link.rsplit('/')[-1]}"
    for link in links:
        fname = get_name(link)
        print(f"downloading {fname}")
        mypool.apply_async(pwget.urlretrieve, args=(f"{web}{link}", fname))

    mypool.close()
    mypool.join()

    # Extract the 'parameters' file, remove the tar and transform it in-place to avoid errors with the gRINN network construction method
    fname = next(get_name(link) for link in links if "prm" in get_name(link))
    with tarfile.open(fname) as tar:
        tar.extractall(path)
    os.remove(fname)

    for line in fileinput.input(f"{path}/parameters", inplace=True):
        if line.strip().startswith("HBOND"):
            line = "HBOND CUTHB 0.5\n"  # gRINN example parameters file has the 'CUTHB 0.5' in addition to 'HBOND' in its line and seems to be necessary
        elif line.strip().startswith("END"):
            line = "END\n"  # gRINN example parameters file has a newline in the end of the file that seems to be necessary
        sys.stdout.write(line)

    return


def process_input(
    GPCR, pdb, protein_sel, pdbf, compdbf, psf, psff, trajs, trajsf, comtrajsf, **kwargs
):
    r"""Process the input structure and trajectory(ies) files

    Only the passed selection is used from the whole input structure file. Its residue
    names are standardized to the usual 3-letter residue codes (
    :func:`AlloViz.AlloViz.trajutils.standardize_resnames`) and, if the structure is a
    GPCR, the generic numbers are retrieved from GPCRdb (
    :func:`AlloViz.AlloViz.trajutils.get_GPCRdb_numbering`). PDB files of the selection
    and its CA atoms are written, and also a PSF file with the selection if applicable.
    The trajectory(ies) of the selection and its CA atoms in xtc format are also
    generated.

    Parameters
    ----------
    GPCR : bool or int
        Use `True` if the structure is a GPCR, or the ID of a GPCRmd database
        dynamics entry if it was used originally to download the files from it.
    pdb : str
        Filename of the PDB structure to read, process and use.
    protein_sel : str
        MDAnalysis atom selection string to select the protein structure from the
        Universe (e.g., in case simulations in biological conditions are used, to avoid
        selecting extra chains, water molecules, ions...).
    pdbf : str
        Complete relative filename of the processed PDB structure to save.
    compdbf : str
        Complete relative filename of the processed PDB structure of the CA atoms to save.
    psf : str
        Filename of the .psf file corresponding to the pdb used to read, process and use.
    psff : str
        Complete relative filename of the processed PSF structure to save.
    trajs : list
        Filename(s) of the MD trajectory (or trajectories) to read and use. File format
        must be recognized by MDAnalysis (e.g., xtc).
    trajsf : list
        Complete relative filename(s) of the processed MD trajectory (or trajectories) to
        save.
    comtrajsf : list
        Complete relative filename(s) of the processed MD trajectory (or trajectories) of
        the CA atoms to save.
    **kwargs
        `special_res` can be passed as an optional kwarg, and it should be a dictionary
        containing a mapping of special residue 3/4-letter code(s) present in the
        structure to the corresponding standard 1-letter code(s).
    """
    whole = mda.Universe(pdb)
    original_resnames = whole.residues.resnames
    protein = whole.select_atoms(protein_sel)

    # Rename all residues in protein_sel to standard names
    standardize_resnames(protein, **kwargs)

    # Retrieve GPCRdb residue generic numbering if it's a GPCR
    if GPCR:
        get_GPCRdb_numbering(protein)

    # If the input files has the same number of atoms and residue names as the processed entities, avoid re-saving pdb and trajectory(ies) files
    if (
        whole.atoms.n_atoms == protein.n_atoms
        and original_resnames == protein.residues.resnames
    ):
        pdbf = pdb
        trajsf = trajs
        psff = psf

    # Make a selection of the CA atoms on the processed protein
    CAs = protein.select_atoms("name CA")

    # Write protein pdb file
    if not os.path.isfile(pdbf):
        protein.write(pdbf)
    if not os.path.isfile(compdbf):
        CAs.write(compdbf)

    # Write protein psf file if it exists
    if psff and not os.path.isfile(psff):
        psf = parmed.load_file(psf)[protein.indices]
        psf.title = psff
        psf.write_psf(psff)

    # Write protein trajectory(ies) file(s): https://stackoverflow.com/a/73043849
    if any(
        [
            not os.path.isfile(f)
            for f in list(trajsf.values()) + list(comtrajsf.values())
        ]
    ):
        whole.load_new(trajs, continuous=False)
        sf = (
            whole.trajectory._start_frames
        )  # for 3 trajs of 2500 frames each looks like: array([   0, 2500, 5000, 7500])
        iterator = enumerate(
            zip(sf[:-1], sf[1:]), 1
        )  # looks like: (1, (0, 2500)), (2, (2500, 5000)), (3, (5000, 7500))

        def write_protein_trajs(start, stop, t, comt):
            if not os.path.isfile(t):
                protein.write(t, frames=whole.trajectory[start:stop])
            if not os.path.isfile(comt):
                # Retrieve the residues' COM trajectory
                COMs = AnalysisFromFunction(lambda atoms: atoms.center_of_mass(compound="residues"),
                                            protein).run(start=start, stop=stop)
                # Create a new universe with the CA atoms as topology and the COMs' trajectory
                COMu = mda.Merge(CAs).load_new(COMs.results['timeseries'], 
                                               format = mda.coordinates.memory.MemoryReader)
                COMu.atoms.write(comt, frames="all")
                

        # For each trajectory asynchronously (with `multiprocess` Pool), write the protein and residues' COM trajectories if the files don't already exist
        mypool = Pool()
        for i, frames in iterator:
            mypool.apply_async(
                write_protein_trajs, args=(*frames, trajsf[i], comtrajsf[i])
            )
        mypool.close()
        mypool.join()

    return pdbf, trajsf, psff



def get_bonded_cys(pdbf):
    r"""Identify disulfide bond-forming cysteines' sulphur atoms

    Returns a list of the residue indices of cysteines whose sulphur atoms form a
    disulfide bond using `ParmEd <https://parmed.github.io/ParmEd/html/index.html>`_.
    These are used downstream by the PyInteraph2_Contacts network construction method to
    delete them from the results, as the program does not detect the bonds by itself and
    it results in contact frequencies of '1' that are uninteresting for allosteric 
    communication.

    Parameters
    ----------
    pdbf : str
        Filename of the already-processed protein PDB structure to use for detection.
    """
    bonded_cys = []
    pdb = parmed.read_pdb(pdbf)
    mask = parmed.amber.mask.AmberMask(pdb, ":CY?@SG")
    for sel in mask.Selected():
        atom = pdb.atoms[sel]
        bonded = [a.name for a in atom.bond_partners]
        if "SG" in bonded:
            bonded_cys.append(atom.residue.idx)

    return bonded_cys
















    # def _get_COM_trajs(self):
    #     compath = f"{self._datadir}/COM_trajs"
    #     if not os.path.isdir(compath): os.makedirs(compath, exist_ok=True)

    #     compdb = f"{compath}/ca.pdb"
    #     if not os.path.isfile(compdb):
    #         print(f"Making trajectories of residue COM for {self.name}")
    #         prot = self.protein.select_atoms("name CA")
    #         prot.write(compdb)
    #     setattr(self, "_compdbf", compdb)

    #     comtrajs = [f"{compath}/{xtc}.xtc" for xtc in self._trajs]
    #     setattr(self, "_comtrajs", {num: traj for num, traj in enumerate(comtrajs, 1)})

    #     for xtc, comtraj in self._comtrajs.items():
    #         if not os.path.isfile(comtraj):
#                 prot = self.protein.atoms
#                 traj = self.u.trajectory.readers[xtc-1] if hasattr(self.u.trajectory, "readers") else self.u.trajectory
#                 # traj = next(traj for traj in self.mdau.trajectory.readers if traj.filename == self._trajs[xtc])
#                 arr = np.empty((prot.n_residues, traj.n_frames, 3))
#                 for ts in traj:
#                     arr[:, ts.frame] = prot.center_of_mass(compound='residues')

#                 cau = mda.Universe(compdb, arr, format=mda.coordinates.memory.MemoryReader, order='afc')
#                 with mda.Writer(comtraj, cau.atoms.n_atoms) as W:
#                     for ts in cau.trajectory:
#                         W.write(cau.atoms)


#     def _translate_ix(self, mapper):
#         # self._translate_ix = lambda mapper: lambda ix: tuple(mapper[_] for _ in ix) if isinstance(ix, tuple) else mapper[ix]
#         return lambda ix: tuple(mapper[_] for _ in ix) if isinstance(ix, tuple) else mapper[ix]

#     def _dihedral_residx(self, end=-1):
#         res_arrays = np.split(self.protein.residues.resindices, np.where(np.diff(self.protein.residues.resnums) != 1)[0]+1)
#         # dihedral_residx = lambda end=-1: [elem for arr in res_arrays for elem in arr[1:end]]
#         return [elem for arr in res_arrays for elem in arr[1:end]]


# def get_mdau(state, protein_sel):
#     pdb = mda.Universe(state._pdbf)#, *state._trajs.values())
#     prot = pdb.select_atoms(protein_sel)
#     mdau = mda.Universe(state._pdbf, *state._trajs.values())

#     protcopy = mda.core.universe.Merge(prot.copy()).atoms
#     # protfile = f"{state._datadir}/prot.pdb"
#     protfile = f"{state._datadir}/prot_{state._pdbf.rsplit('/', 1)[-1]}"
#     setattr(state, "_protpdb", protfile)

#     try:
#         from Bio.SeqUtils import seq1, seq3
#         for res in protcopy.residues:
#             res.resname = seq3(seq1(res.resname, custom_map = state._res_dict)).upper()
#     except:
#         print(f"""
# resname(s) {[res for res in set(protcopy.residues.resnames) if res not in state._res_dict]} could not be translated to a standard name.
# Please provide a mapping from 3/4-letter code to 1-letter code as a dictionary with keyword argument 'special_res' upon initialization.
#               """)
#         raise


#     if state.GPCR:#hasattr(state, "_gpcrmdid"):
#         # prot_numsf = f"{state._datadir}/{state._pdbf.replace(f'{state._path}/', 'gpcrdb_')}"
#         prot_numsf = state._protpdb

#         if not os.path.isfile(prot_numsf):
#             print(f"retrieving {prot_numsf}")
#             with io.StringIO() as protf:
#                 with mda.lib.util.NamedStream(protf, "prot.pdb") as f:
#                     protcopy.write(f, file_format="pdb")
#                 response = requests.post('https://gpcrdb.org/services/structure/assign_generic_numbers',
#                                          files = {'pdb_file': protf.getvalue()})
#                 with open(prot_numsf, "w") as prot_nums:
#                     prot_nums.write(response.text)

#         # print "GPCRdb script labeling generic numbers on the annotated pdb structure\nKey bindings for labels:\nF1 - show generic numbers\nF2 - show Ballesteros-Weinstein numbers\nF3 - clear labels"
#         # label n. CA or n. N
#         # cmd.set_key('F1', 'label n. CA & (b >-8.1 and  b <8.1), str("%1.2f" %b).replace(".","x") if b > 0 else str("%1.3f" %(-b + 0.001)).replace(".", "x")')
#         # cmd.set_key('F2', 'label n. N & (b > 0 and  b <8.1), "%1.2f" %b')
#         # cmd.set_key('F3', 'label n. CA or n. N')
#         nums = mda.Universe(prot_numsf).select_atoms(f"{protein_sel} and name CA").tempfactors
#         GPCRdb_nums = np.array([round(bfactor, 2) if bfactor>1 else round(-bfactor+0.001, 3) if bfactor<0 else 0 for bfactor in nums])

#         mdau.select_atoms(f"{protein_sel} and name CA").tempfactors = GPCRdb_nums#.round(2)
#         pdb.select_atoms(f"{protein_sel} and name CA").tempfactors = GPCRdb_nums
#         prot.select_atoms(f"{protein_sel} and name CA").tempfactors = GPCRdb_nums

#     else:
#         protcopy.write(state._protpdb)


#     return pdb, prot, mdau


# def get_dihedral_residx(prot):
#     res_arrays = np.split(prot.residues.resindices, np.where(np.diff(prot.residues.resnums) != 1)[0]+1)
#     dihedral_residx = lambda end=-1: [elem for arr in res_arrays for elem in arr[1:end]]

#     return dihedral_residx


# def get_res_dict(state, **kwargs):
#     from Bio.SeqUtils import IUPACData
#     from MDAnalysis.lib.util import inverse_aa_codes

#     res_d = {}
#     res_d.update(IUPACData.protein_letters_3to1_extended)
#     res_d.update(inverse_aa_codes)
#     if "special_res" in kwargs and isinstance(kwargs["special_res"], dict):
#         res_d.update(kwargs["special_res"])

#     return res_d


#######3 change state to Protein and use Protein's "u"; don't use protein sel


# def get_bonded_cys(state, pdb):
#     import parmed

#     indices = []

#     pdb = parmed.read_pdb(pdb)
#     mask = parmed.amber.mask.AmberMask(pdb, ':CY?@SG')
#     for sel in mask.Selected():
#         atom = pdb.atoms[sel]
#         bonded = [a.name for a in atom.bond_partners]
#         if "SG" in bonded:
#             indices.append(atom.residue.idx)

#     return indices
