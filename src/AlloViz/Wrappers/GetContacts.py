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
"_getcontacts_contacts": ".Packages.getcontacts.get_dynamic_contacts",
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
        ctcs = f"/home/{os.environ.get('USER')}/{self._d['name']}/{xtc}.tsv"
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
        df.index = df.index.map(lambda idx: tuple(sorted(idx, key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(self._rawpq(xtc))