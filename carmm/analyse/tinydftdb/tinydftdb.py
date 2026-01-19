from tinydb import TinyDB, Query, table
import os

class TinyDFTDB():

    def __init__(self, dbdir="./", PID=None, save_structure=True, struct_dir=None, human_readable_hash=True):

        self.db_pid = PID
        self.dbdir = dbdir
        self.human_readable_hash = human_readable_hash
        self.save_structure = save_structure

        if struct_dir is None:
            self.struct_dir = f"{self.dbdir}/{self.db_pid}_STRUCTURES"
        else:
            self.struct_dir = struct_dir

        if any(self.db_pid) is None:
            raise Exception("Project ID must be specified")

        self.db = TinyDB(f'{self.dbdir}/{self.db_pid}.json')

    def insert_calc_start_dataline(self, platform, EID, DID, CID, **kwargs):
        
        self.table = self.db.table(EID)
        self.table.document_id_class = str

        meta_data = self.generate_metadata(platform, EID, DID, CID, **kwargs)
        self.init_dict = meta_data | {"InputParameters": self.generate_calulcation_inputs(**kwargs)} | {"TimingData" : self.timing_data}

        self.table = self.db.table(EID)

        self.init_dict["CalculationStatus"] = "In Progress"

        self.insert_dataline(self.init_dict)

    def insert_calc_end_dataline(self, calc_success=True, **kwargs):
        import datetime
        import os
        from ase.io import write

        if calc_success:
            self.init_dict["CalculationStatus"] = "Success"            
        else:
            self.init_dict["CalculationStatus"] = "Failed"

        self.end_datetime = datetime.datetime.now()
        self.init_dict["TimingData"]["EndDate"] = self.end_datetime.strftime("%d%m%y")
        self.init_dict["TimingData"]["EndTime"] = self.end_datetime.strftime("%X")
        delta_calc = self.end_datetime - self.start_datetime
        self.init_dict["TimingData"]["ElapsedTime"] = delta_calc.total_seconds()

        atom_keys = []
        for key, val in kwargs.items():
            if type(val).__name__ == "Atoms":
                atom_keys.append(key)
                if (self.save_structure):
                    os.makedirs(self.struct_dir, exist_ok=True)
                    write(f"{self.struct_dir}/FINAL_{self.uid_hash}.traj", val)

        for key in atom_keys:
            _ = kwargs.pop(key)

        if len(kwargs) > 0:
            self.final_dict = self.init_dict | {"CalculationOutput": kwargs}
        else:
            self.final_dict = self.init_dict
            
        self.insert_dataline(self.final_dict)

    def insert_dataline(self, data):
        self.table.upsert(table.Document(data, doc_id=self.uid_hash))
        
    def generate_uid(self, EID, DID, CID):
        from hashlib import sha256

        if self.human_readable_hash:
            uid_hash = "_".join([self.db_pid, EID, DID, CID])
        else:
            uid_str = "".join([self.db_pid, EID, DID, CID]).encode('UTF-8')
            uid_hash = sha256(uid_str).hexdigest()

        return uid_hash

    def generate_metadata(self, platform, EID, DID, CID, **kwargs):
        import datetime

        data_entry = {}

        # Define metadata
        metadata = {}
        metadata["PID"] = self.db_pid
        self.EID = EID
        self.DID = DID
        self.CID = CID
        metadata["EID"] = self.EID
        metadata["DID"] = self.DID
        metadata["CID"] = self.CID

        self.start_datetime = datetime.datetime.now()
        self.timing_data = {}
        self.timing_data["StartDate"] = self.start_datetime.strftime("%d%m%y")
        self.timing_data["StartTime"] = self.start_datetime.strftime("%X")

        data_entry["Metadata"]  = metadata

        # HPC metadata - we're going to make a wild assumption that
        # anything that isn't a desktop is a SLURM based HPC
        hpc_data_names = ["Platform", "SLURM_JOB_ID", "SLURM_SUBMIT_DIR", "SLURM_JOB_NODELIST", "SLURM_NTASKS", "SLURM_NTASKS_PER_NODE"]
        hpc_data_entries = {}

        hpc_data_entries["Platform"] = platform
        for name in hpc_data_names[1:]:
            if name in os.environ:
                hpc_data_entries[name] = os.environ[name]

        data_entry["HPCMetadata"] = hpc_data_entries

        self.uid_hash = self.generate_uid(EID, DID, CID)

        return data_entry

    def generate_calulcation_inputs(self, **kwargs):
        import os
        from ase.io import write

        from ase.calculators.aims import Aims, AimsProfile
        from ase.calculators.genericfileio import GenericFileIOCalculator

        input_param_dict = {}

        for key, val in kwargs.items():

            if hasattr(val, "__class__"):
                if issubclass(val.__class__, GenericFileIOCalculator):
                    input_param_dict[key] = val.parameters
                    input_param_dict[key] = input_param_dict[key] | {"calc_name":type(val).__name__}

                if type(val).__name__ == "Atoms":
                    at_dict = {}
                    at_dict = at_dict | {"Formula" : val.get_chemical_formula()}
                    at_dict = at_dict | {"NAtoms" : len(val)}
                    input_param_dict[key] = at_dict
                    if self.save_structure:
                        os.makedirs(self.struct_dir, exist_ok=True)
                        write(f"{self.struct_dir}/INIT_{self.uid_hash}.traj", val)
            else:
                input_param_dict[key] = val

        return input_param_dict

        
if __name__ == "__main__":

    from ase.calculators.aims import Aims, AimsProfile
    from ase.calculators.genericfileio import GenericFileIOCalculator
    from ase import Atoms

    atoms=Atoms("H2O")

    calc = Aims(profile=AimsProfile(command="NoCommand"),
                xc='LDA',
                output=['dipole'],
                sc_accuracy_etot=1e-2,
                sc_accuracy_eev=1e-1,
                sc_accuracy_rho=1e-2,
                sc_accuracy_forces=1e-1,)

    myfirstdb = TinyDFTDB(PID="BMARK")

    myfirstdb.insert_calc_start_dataline(platform="desktop",
                                         EID="test_exp",
                                         DID="test_conv",
                                         CID="test_param",
                                         atoms=atoms,
                                         calculator=calc)

    myfirstdb.insert_calc_end_dataline(calc_success=True,
                                       total_energy=5.5)
    
