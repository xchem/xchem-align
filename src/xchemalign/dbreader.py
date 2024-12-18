# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sqlite3
import pandas as pd

# RefinementOutcome values:
# 1 - Analysis Pending
# 2 - PANDDA model
# 3 - In Refinement
# 4 - CompChem ready
# 5 - Deposition ready
# 6 - Deposited
# 7 - Analysed & Rejected


def read_dbmeta(dbfile):
    """
    Read the XCE metadata from the mainTable table of a sqllite db
    :param dbfile: The sqlite db file
    :return: Pandas dataframe with the contends of the table
    """
    # Create your connection.
    cnx = sqlite3.connect(dbfile)
    df = pd.read_sql_query(
        """SELECT ID, CompoundSMILES, CompoundCode, CrystalName, ispybStatus,
                            RefinementCIF, RefinementCIFStatus, RefinementBoundConformation, RefinementMTZ_latest,
                            RefinementDate, RefinementOutcome, LastUpdated
                            FROM mainTable WHERE RefinementOutcome IS NOT NULL""",
        cnx,
    )

    df["LastUpdatedDate"] = pd.to_datetime(df["LastUpdated"], infer_datetime_format=True)
    return df


def filter_dbmeta(dbfile, reference_datasets):
    df1 = read_dbmeta(dbfile)
    df2 = df1[
        (
            df1.CrystalName.isin(reference_datasets)
            | df1.RefinementOutcome.str.startswith("4")
            | df1.RefinementOutcome.str.startswith("5")
            | df1.RefinementOutcome.str.startswith("6")
            | df1.RefinementOutcome.str.startswith("7")
        )
    ]
    return df2
