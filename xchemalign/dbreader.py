import sqlite3
import pandas as pd


def read_dbmeta(dbfile):
    """
    Read the XCE metadata from the mainTable table of a sqllite db
    :param dbfile: The sqllite db file
    :return: Pandas dataframe with the contends of the table
    """
    # Create your connection.
    cnx = sqlite3.connect(dbfile)
    df = pd.read_sql_query('''SELECT ID, CompoundSMILES, CompoundCode, CrystalName, ispybStatus, 
                            RefinementCIF, RefinementCIFStatus, RefinementPDB_latest, RefinementMTZ_latest,
                            RefinementDate, RefinementOutcome
                            FROM mainTable WHERE RefinementOutcome IS NOT NULL''', cnx)
    return df


def filter_dbmeta(dbfile):
    df1 = read_dbmeta(dbfile)
    df2 = df1[(
               df1.RefinementOutcome.str.startswith('3') |
               df1.RefinementOutcome.str.startswith('5') |
               df1.RefinementOutcome.str.startswith('6'))]
    return df2


def main():

    df = filter_dbmeta('data/dls/labxchem/data/lb18145/lb18145-216/processing/database/soakDBDataFile.sqlite')
    print(df.RefinementPDB_latest)


if __name__ == "__main__":
    main()