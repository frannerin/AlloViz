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
    def _computation(self, xtc):
        """"""
        path = self._path
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        if not os.path.isfile(freqs):
            _getcontacts_contacts.main(f"--topology {self._pdbf} --trajectory {self._trajs[xtc]} --output {ctcs} --itypes all --cores {self.taskcpus}".split())
            _getcontacts_freqs.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc
        
        
    def _save_pq(self, args):
        freqs, xtc = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        df.index = df.index.map(lambda idx: tuple(sorted([res.split(":", 1)[-1] for res in idx], key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(self._rawpq(xtc))