import re
import os
import duckdb
class Feature:

    def __init__(self, parquet_path: str):
        if os.path.exists(parquet_path):
            self.parquet_db = duckdb.connect()
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(parquet_path))
        else:
            raise FileNotFoundError(f'the file {parquet_path} does not exist.')

    def get_report_from_database(self, samples: list):
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param runs: A list of ms_runs
        :return: The report
        """
        database = self.parquet_db.sql(
                        """
            select * from parquet_db
            where sample_accession IN {}
            """.format(tuple(samples))
        ) 
        report = database.df()
        return report

    def iter_samples(self,file_num:int=20):
        """
        :params file_num: The number of files being processed at the same time(default 10)
        :yield: _description_
        """
        samples = self.get_unique_samples()
        ref_list =  [samples[i:i+file_num] for i in range(0,len(samples), file_num)]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs)
            yield refs,batch_df

    def get_unique_samples(self):
        """
        return: A list of deduplicated peptides.
        """
        unique_peps = self.parquet_db.sql(f"SELECT DISTINCT sample_accession FROM parquet_db").df()

        return unique_peps['sample_accession'].tolist()