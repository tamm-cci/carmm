from tinydb import TinyDB, Query, table
import os

class TinyDFTDB():

    def __init__(self, dbdir="./", PID=None):

        self.db_pid = PID
        self.dbdir = dbdir

        if any(self.db_pid) is None:
            raise Exception("Project ID must be specified")

        self.db = TinyDB(f'{self.dbdir}/{self.db_pid}.json')

    def insert_calc_start_dataline(self, platform, EID, DID, CID, **kwargs):
        self.table = self.db.table(EID)
        self.table.document_id_class = str

        meta_data = self.generate_metadata(platform, EID, DID, CID, **kwargs)
        self.init_dict = meta_data | {"input_parameters": self.generate_calulcation_inputs(**kwargs)}

        self.table = self.db.table(EID)

        self.init_dict["CalculationStatus"] = "In Progress"

        self.insert_dataline(self.init_dict)

    def insert_calc_end_dataline(self, calc_success=True, **kwargs):

        if calc_success:
            self.init_dict["CalculationStatus"] = "Success"
        else:
            self.init_dict["CalculationStatus"] = "Failed"

        if len(kwargs) > 0:
            self.final_dict = self.init_dict | {"CalculationOutput": kwargs}
        else:
            self.final_dict = self.init_dict
            
        self.insert_dataline(self.final_dict)

    def insert_dataline(self, data):
        self.table.upsert(table.Document(data, doc_id=self.uid_hash))
        
    def generate_uid(self, EID, DID, CID):
        from hashlib import sha256

        uid_str = "".join([self.db_pid, EID, DID, CID]).encode('UTF-8')
        uid_hash = sha256(uid_str).hexdigest()

        return uid_hash

    def generate_metadata(self, platform, EID, DID, CID, **kwargs):
        import datetime

        data_entry = {}
        # Define metadata
        metadata = {}
        metadata["PID"] = self.db_pid
        metadata["EID"] = EID
        metadata["DID"] = DID
        metadata["CID"] = CID

        curr_datetime = datetime.datetime.now()

        metadata["Date"] = curr_datetime.strftime("%d%m%y")
        metadata["StartTime"] = curr_datetime.strftime("%X")

        data_entry["metadata"]  = metadata

        # HPC metadata - we're going to make a wild assumption that
        # anything that isn't a desktop is a SLURM based HPC
        if platform != "desktop":
            hpc_data_names = ["Platform", "JobStatus", "JOB_ID", "SUBMIT_DIR", "NODELIST", "NTASKS", "NTASKSPERNODE"]
            hpc_data_entries = {}
            hpc_data_entries[hpc_data_names[0]] = platform
            hpc_data_entries[hpc_data_names[1]] = os.environ["SLURM_JOB_ID"]
            hpc_data_entries[hpc_data_names[2]] = os.environ["SLURM_SUBMIT_DIR"]
            hpc_data_entries[hpc_data_names[3]] = os.environ["SLURM_JOB_NODELIST"]
            hpc_data_entries[hpc_data_names[4]] = os.environ["SLURM_NTASKS"]
            hpc_data_entries[hpc_data_names[5]] = os.environ["SLURM_NTASKS_PER_NODE"]

            data_entry["hpc_metadata"] = hpc_data_entries

        self.uid_hash = self.generate_uid(EID, DID, CID)

        return data_entry

    def generate_calulcation_inputs(self, **kwargs):

        print(kwargs.values())

        input_param_dict = {}

        for key, val in kwargs.items():

            if hasattr(val, "__class__"):
                if issubclass(val.__class__, GenericFileIOCalculator):
                    input_param_dict[key] = val.parameters
                    input_param_dict[key] = input_param_dict[key] | {"calc_name":type(val).__name__}

                if type(val).__name__ == "Atoms":
                    at_dict = {}
                    at_dict = at_dict | {"formula" : val.get_chemical_formula()}
                    at_dict = at_dict | {"natoms" : len(val)}
                    input_param_dict[key] = at_dict

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
    
