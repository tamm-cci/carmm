from tinydb import Query, TinyDB

class TinyDFTDB_Explorer():

    def __init__(self, db_dir, PID):
        self.db = TinyDB(f'{db_dir}/{PID}.json')

    @property
    def table_names(self):
        return sorted(list(self.db.tables()))

    def load_table(self, uid):
        if type(uid) is int:
            uid = self.table_names[uid]

        table = self.db.table(uid)
        table.document_id_class = str
        return table

    def print_key_tree(self, table, nid):
        table_dict = table.all()[nid]

        def unpack_layer(curr_dict, print_str, idx):
            try:
                if type(curr_dict[idx]) is dict:
                    print_str=f"|    {print_str}"
                    for idx1 in curr_dict[idx]:
                        print(f"{print_str} {idx1}")
                        unpack_layer(curr_dict[idx], print_str, idx1)
            except:
                return None

        print(f"PRINTING KEY TREE UP TO THREE LAYERS")
        print_str = "|--"
        for idx in table_dict:
            print(f"{print_str} {idx}")            
            unpack_layer(table_dict, print_str, idx)

    def load_columns_to_pandas(self, table, ind_vars, *args):
        from glom import glom
        import pandas as pd

        if type(ind_vars) == list:
            columns = ind_vars + list(args)
        else:
            columns = [ind_vars] + list(args)

        short_columns = []
        for name in columns:
            short_columns.append(name.split(".")[-1])

        data_dict = {}
        for idx, col in enumerate(columns):
            data_dict[short_columns[idx]] = []
            for row in table.all():
                try:
                    data_dict[short_columns[idx]].append(glom(row, col))
                except:
                    data_dict[short_columns[idx]].append("NaN")
                    
        return pd.DataFrame(data_dict)

    def print_calc_status(self):

        for i in range(len(db_exp.table_names)):

            print(f"|------ EXPERIMENT:{db_exp.table_names[i]} STATUS -------|")
            table = db_exp.load_table(i)
            df = db_exp.load_columns_to_pandas(table, "Metadata.DID", "Metadata.CID",
                                           "TimingData.StartTime", "CalculationStatus")

            df = df.sort_values(by=["DID","CID"])

            print(df)
            print("\n")

if __name__ == "__main__":
    #
    #  A funky little script for outputting all calculations and
    #  their progress.
    #
    #  load_columns_to_pandas is a little smarter and loads data
    #  into a pandas dataframe for post-processing
    #
    
    from glom import glom
    
    db_exp = TinyDFTDB_Explorer("./ProjectData", "Benchmarking")

    print(db_exp.table_names)
    table = db_exp.load_table(0)

    db_exp.print_calc_status()
    
