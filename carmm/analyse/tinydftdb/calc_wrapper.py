from ase.calculators.aims import Aims, AimsProfile
from ase.calculators.emt import EMT
from ase.calculators.genericfileio import GenericFileIOCalculator
from carmm.analyse.tinydftdb.tinydftdb import TinyDFTDB
from ase import Atoms

class TinyDFTDBCalcWrapper():

    def __init__(self, PID, db_dir="./"):
        from carmm.analyse.tinydftdb.tinydftdb import TinyDFTDB

        self.PID = PID
        self.db_dir = db_dir

        self.dftdb = TinyDFTDB(PID=PID, save_structure=True)


    def init_dataline(self, EID, DID, CID, platform, **kwargs):

        self.dftdb.insert_calc_start_dataline(platform=platform,
                                             EID=EID,
                                             DID=DID,
                                             CID=CID,
                                             **kwargs)


    def set_calc_function(self, calc_function):
        self.calc_function = calc_function


    def run_calculation(self, *args, **kwargs):
        import traceback
        import tracemalloc
        import sys

        try:
            tracemalloc.start()
            output_dict = self.calc_function(*args, **kwargs)
        except Exception as eee:
            curr_size, curr_peak = tracemalloc.get_traced_memory()
            self.dftdb.insert_calc_end_dataline(calc_success=False, PeakMemoryMB=curr_peak/(1024*1024))
            print(f"Something has happened while running the calculation: {eee}\nAborting...")
            traceback.print_tb(eee.__traceback__, limit=5, file=sys.stdout)
        else:
            curr_size, curr_peak = tracemalloc.get_traced_memory()
            output_dict = {"PeakMemoryMB": curr_peak/(1024*1024)} | output_dict
            self.dftdb.insert_calc_end_dataline(calc_success=True, **output_dict)

if __name__ == "__main__":
    atoms1=Atoms("H2O")
    calc = EMT()
    atoms1.calc = calc

    atoms2=Atoms("Pu")
    calc = EMT()
    atoms2.calc = calc

    
    def run_function(atoms):
        tote = atoms.get_potential_energy()
        data_dict = {"TotalEnergy": tote}

        return data_dict

    db_wrapper = TinyDFTDBCalcWrapper(PID="Benchmarking")
    db_wrapper.set_calc_function(run_function)

    static_metadata = {"platform":"Desktop", "EID":"BasisConvergence", "DID":"Light"}

    db_wrapper.init_dataline(**static_metadata, CID="H2O", atoms=atoms1, calculator=calc)
    db_wrapper.run_calculation(atoms1)

    db_wrapper.init_dataline(**static_metadata, CID="H2", atoms=atoms2, calculator=calc)
    db_wrapper.run_calculation(atoms2)
