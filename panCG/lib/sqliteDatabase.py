import sqlite3


class SQLiteDatabase:
    def __init__(self, db_name: str, index_column: str, data_columns: list):
        """
        Initialize the database connection and specify the index column and data column
        Args:
            db_name: SQLite database file name
            index_column: The column name to use as the index
            data_columns: List of column names of data columns
        """
        self.db_name = db_name
        self.index_column = index_column
        self.data_columns = data_columns

    @staticmethod
    def _get_sql_type(column):
        """
        Automatically select the data type supported by SQLite based on the data type of the column
        Args:
            column: Pandas columns
        Returns:
            SQLite Data Types
        """
        if column.dtype == 'int64':
            return 'INTEGER'
        elif column.dtype == 'float64':
            return 'REAL'
        elif column.dtype == 'object':
            return 'TEXT'
        else:
            return 'TEXT'

    def insert_data(self, dfs: list):
        """
        Writing multiple DataFrames to a SQLite database
        Args:
            dfs: List of DataFrames
        Returns:
        """
        conn = sqlite3.connect(self.db_name)
        cursor = conn.cursor()
        for df in dfs:
            column_defs = ', '.join([f'{col} {self._get_sql_type(df[col])}' for col in df.columns])
            create_table_sql = f"CREATE TABLE IF NOT EXISTS data ({column_defs})"
            cursor.execute(create_table_sql)
            df.to_sql('data', conn, if_exists='append', index=False)  # Writing the DataFrame to a Database
        conn.commit()  # Committing a transaction
        conn.close()   # Close the connection

    def create_index(self):
        """
        Create Index
        Returns:
        """
        conn = sqlite3.connect(self.db_name)
        cursor = conn.cursor()
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{self.index_column} ON data({self.index_column})")
        conn.commit()
        conn.close()

    def query_data(self, ref_to_find):
        """
        Query the data corresponding to the specified ref
        Args:
            ref_to_find: Query the ref value
        Returns:
            Query results
        """
        conn = sqlite3.connect(self.db_name)
        # conn.execute('PRAGMA journal_mode=WAL;')  # Enable WAL mode
        cursor = conn.cursor()
        columns_str = ', '.join(self.data_columns)  # Get all data columns
        cursor.execute(f"SELECT {columns_str} FROM data WHERE {self.index_column} = ?", (ref_to_find,))
        result = cursor.fetchall()
        conn.close()
        return result

    def commit_and_close(self):
        """
        Commit the transaction and close the database connection
        """
        conn = sqlite3.connect(self.db_name)
        conn.commit()
        conn.close()

