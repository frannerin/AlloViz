"""GetContacts wrapper

It calculates contact frequencies.

"""

import os, time, sys

from pkgutil import get_loader

import pandas

from .Base import lazy_import, Multicore

from ..AlloViz.utils import get_pool


sys.path.append(os.path.dirname(get_loader("AlloViz").path) + "/Packages/getcontacts")

imports = {
# "_getcontacts_contacts": ".Packages.getcontacts.get_dynamic_contacts",
"_getcontacts_freqs": ".Packages.getcontacts.get_contact_frequencies",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    


class GetContacts(Multicore):
    """GetContacts' contact frequencies
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        if "GetContacts_threshold" in d:
            new.GetContacts_threshold = d["GetContacts_threshold"]
        return new     
    
    def __init__(self, *args):
        super().__init__(*args)
        
        # Filter dataset according to GetContacts_threshold optional kwarg
        if hasattr(self, "GetContacts_threshold"):
            # Define the list of .pq files that we expect are going to be saved (or be retrieved) and a function to check which of them already exist
            pqs = [self._rawpq(xtc) for xtc in self._trajs]
            no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]

            # Function to wait for the calculations to finish in the background; returns the .pq files to be read and added as attributes when they do
            def wait_calculate(pqs):
                while any(no_exist(pqs)):
                    time.sleep(5)
                return pqs
            
            # Wait asynchronously for analysis to end and then add the filter data
            filter_raw = lambda _: setattr(self, "raw", self.raw[self.raw["weight"] >= self.GetContacts_threshold])
            get_pool().apply_async(wait_calculate,
                                   args=(pqs,),
                                   callback=filter_raw)
            
            
    def _computation(self, xtc):
        """"""
        path = self._path
        #ctcs = f"{path}/{xtc}.tsv"
        ctcs = f"/users/gpcr/{os.environ.get('USER')}/{self._d['name']}/{xtc}.tsv"
        print(ctcs, os.path.isfile(ctcs))
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        if not os.path.isfile(freqs):# or ow:
            #_getcontacts_contacts.main(f"--topology {self._pdbf} --trajectory {self._trajs[xtc]} --output {ctcs} --itypes all --cores {self.taskcpus}".split())
            _getcontacts_freqs.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc
        
        
    def _save_pq(self, args):
        freqs, xtc = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        # df.index = df.index.map(lambda idx: tuple(sorted([res.split(":", 1)[-1] for res in idx], key = lambda res: int(res.split(":")[-1]))))
        
        # Filter out rows that contain contacts with the ligand; ligand resid won't be in the Protein's "protein" attribute
        df = df[[all(int(res.split(":")[-1]) in self._d["protein"].residues.resids for res in ix) for ix in df.index]]
        
        # Added to process residue 3-letter codes to change them from the ones in GPCRmd to the standard form that the files used by AlloViz have
        from Bio.SeqUtils import seq1, seq3
        process = lambda name: seq3(seq1(name, custom_map=self._d["_standard_resdict"])).upper()
        df.index = df.index.map(lambda idx: tuple(":".join([chain, process(name), num]) for res in idx for chain, name, num in [res.split(":")]))
        
        df.index = df.index.map(lambda idx: tuple(sorted(idx, key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(self._rawpq(xtc))