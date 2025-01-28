import os

from threading import Thread
from queue import Queue, Empty
from typing import Any

import pandas as pd
import pyarrow as pa

from pyarrow import parquet as pq


class WriteCSVThread(Thread):
    path: str

    write_options: dict[str, Any]

    _queue: Queue
    _wrote_header: bool

    def __init__(
        self, path: str, daemon: bool = True, write_options: dict | None = None, **kwargs
    ):
        super().__init__(daemon=daemon)
        if write_options is None:
            write_options = {}

        path, _ext = os.path.splitext(path)
        path += ".csv"

        self.path = path
        self.write_options = write_options | kwargs
        self._wrote_header = False
        self._queue = Queue()

    def write(self, table: pd.DataFrame):
        self._queue.put(table)

    def close(self):
        self._queue.put(None)
        self.join()

    def run(self):
        while True:
            try:
                table: pd.DataFrame | None = self._queue.get(True)
            except Empty:
                continue

            if table is None:
                break

            table.to_csv(
                self.path,
                header=not self._wrote_header,
                mode="a+" if self._wrote_header else "w",
                index=False,
                **self.write_options
            )
            self._wrote_header = True


class WriteParquetThread(Thread):
    path: str
    metadata: dict[str, Any]

    _queue: Queue

    _schema: pa.Schema | None
    _writer: pq.ParquetWriter | None

    def __init__(self, path: str, daemon: bool = True, metadata: dict | None=None, **kwargs):
        super().__init__(daemon=daemon)

        if metadata is None:
            metadata = {}
        path, _ext = os.path.splitext(path)
        path += '.parquet'

        self.path = path
        self.metadata = metadata | kwargs
        self._queue = Queue()
        self._writer = None
        self._schema = None

    def write(self, table: pd.DataFrame):
        self._queue.put(table)

    def close(self):
        self._queue.put(None)
        self.join()

    def run(self):
        while True:
            try:
                table: pd.DataFrame | None = self._queue.get(True)
            except Empty:
                continue

            if table is None:
                break

            if self._schema is None:
                self._schema = pa.Schema.from_pandas(table, preserve_index=False)
                self._writer = pq.ParquetWriter(self.path, schema=self._schema)

            arrow_table = pa.Table.from_pandas(table, preserve_index=False)
            self._writer.write_table(arrow_table)

        self._writer.add_key_value_metadata(self.metadata)
        self._writer.close()
