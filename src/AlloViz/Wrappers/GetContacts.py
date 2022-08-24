"""GetContacts wrapper

It calculates contact frequencies.

"""

import os

import pandas

from .Base import lazy_import, Multicore

imports = {
"_getcontacts_contacts": ".Packages.getcontacts.get_dynamic_contacts",
"_getcontacts_freqs": ".Packages.getcontacts.get_contact_frequencies",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    


class GetContacts(Multicore):
    """GetContacts' contact frequencies
    """
            
    def _computation(self, xtc):#pdb, traj, xtc, pq, ctcs, freqs, taskcpus):
        """"""
        path = self._path
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        if not os.path.isfile(freqs):# or ow:
            _getcontacts_contacts.main(f"--topology {self._pdbf} --trajectory {self._trajs[xtc]} --output {ctcs} --itypes all --cores {self.taskcpus}".split())
            _getcontacts_freqs.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc
        
        
    def _save_pq(self, args):
        freqs, xtc = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        df.index = df.index.map(lambda idx: tuple(sorted([res.split(":", 1)[-1] for res in idx], key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(self._rawpq(xtc))
        
    
    
    def __init__(self, *args):
        super().__init__(*args)
        
        # Filter dataset according to GetContacts_threshold optional kwarg
        if "GetContacts_threshold" in self._d:
            self.filter_contacts(self._d["GetContacts_threshold"])
    
    @staticmethod
    def _filter_raw(raw, GetContacts_threshold):
        """"""
        return raw[raw["weight"] >= GetContacts_threshold]
    
    def filter_contacts(self, GetContacts_threshold:float):
        r"""Filter contacts below a frequency threshold.
        
        Method.

        Parameters
        ----------
        GetContacts_threshold : float
            Value of the minimum contact frequency (between 0 and 1) threshold, which
            will be used to filter out contacts with a frequency (average) lower than it.

        Returns
        -------
        None
                                    
        Notes
        -----
        Method returns nothing, but GetContacts raw data attribute is filtered according
        to the threshold. Make sure to delete previous analysis attributes and files.
        Analysis results are determined by the input data (e.g., filtered or unfiltered
        network), so they can't simply be filtered according to this new criteria and
        the analysis must be done again.
        """
        print("Make sure to delete/have deleted all previous analysis attributes and files.")
        # message = True
        for filterby in utils.filterbysl:
            if hasattr(self, filterby.capitalize()):
                # if message:
                #     print("It looks like you had already analyzed the network before this contacts filtering, please make sure you delete all analysis files so that the analysis is done again with the filtered network.")
                # message = False
                
                delattr(self, filterby.capitalize())
        
        
        self.raw = self._filter_raw(self.raw, GetContacts_threshold)
        
        