from Pkgs import *

sys.path.append('getcontacts')
from getcontacts import get_dynamic_contacts, get_contact_frequencies


class Getcontacts(Multicorepkg):
    def __init__(self, state):
        # try:
        #print(sys.path)
        #sys.path.append('getcontacts')
        #from getcontacts import get_dynamic_contacts, get_contact_frequencies
        # except:
        #     print("Can't import getcontacts module")
        
        super().__init__(state)
        
        
        
        
        
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        
        path = self._path("raw")
        ctcs = f"{path}/{xtc}.tsv"
        freqs = f"{path}/{xtc}_freqs.tsv"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, xtc, pq, ctcs, freqs, self.taskcpus),
                         callback=self._save_pq)
        
        for _ in range(self.taskcpus-1): pool.apply_async(self._calculate_empty, args=(pq,))
        
        
    def _computation(self, pdb, traj, xtc, pq, ctcs, freqs, taskcpus):
        print("computing")
        if not os.path.isfile(freqs):# or ow:
            get_dynamic_contacts.main(f"--topology {pdb} --trajectory {traj} --output {ctcs} --itypes all --cores {taskcpus}".split())
            get_contact_frequencies.main(f"--input_files {ctcs} --output_file {freqs}".split())
        return freqs, xtc, pq
        
        
    def _save_pq(self, args):
        freqs, xtc, pq = args
        
        df = pandas.read_csv(freqs, sep="\t", skiprows=2,
                             index_col = (0, 1), names = [f"{xtc}"])
        df.index = df.index.map(lambda idx: tuple(sorted(idx, key = lambda res: int(res.split(":")[-1]))))
        df.to_parquet(pq)
