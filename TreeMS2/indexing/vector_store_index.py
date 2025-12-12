import math
import os
import time
import tracemalloc
from typing import Tuple, Optional, Iterator

import faiss
import numpy as np
import psutil
from faiss.contrib import clustering
from tqdm import tqdm

from TreeMS2.config.env_variables import TREEMS2_NUM_CPUS, TREEMS2_MEM_PER_CPU
from TreeMS2.config.logger_config import get_logger, format_execution_time
from TreeMS2.ingestion.storage.vector_store import VectorStore

logger = get_logger(__name__)


class IndexConstructionError(Exception):
    pass


class VectorStoreIndex:
    def __init__(self, vector_store: VectorStore):
        """
        Index for fast ms/ms ingestion similarity search.
        """
        self.vector_store = vector_store
        self.vector_count = vector_store.count_vectors()

        factory_string, nlist = self._create_factory_string(
            vector_count=self.vector_count, vector_dim=vector_store.vector_dim
        )
        self.index = faiss.index_factory(
            vector_store.vector_dim, factory_string, faiss.METRIC_INNER_PRODUCT
        )
        self.factory_string = factory_string
        self.nlist = nlist

    @staticmethod
    def _estimate_memory_budget():
        system_total_ram = psutil.virtual_memory().total
        current_process_ram_usage = psutil.Process().memory_info().rss

        # calculate the amount of RAM that will be reserved for the OS
        giga_byte = 1024**3
        system_total_ram_gb = system_total_ram // giga_byte
        ram_reserved_for_os = 2 * giga_byte  # initial RAM reserved
        if system_total_ram > 4 * giga_byte:
            ram_reserved_for_os += (
                (min(system_total_ram_gb, 16) - 4) // 4
            ) * giga_byte  # +1GB per 4GB (up to 16GB)
        if system_total_ram > 16 * giga_byte:
            ram_reserved_for_os += (
                (system_total_ram_gb - 16) // 8
            ) * giga_byte  # +1GB per 8GB above 16GB

        memory_budget = max(
            system_total_ram - current_process_ram_usage - ram_reserved_for_os, 0
        )
        memory_budget = memory_budget / 2.0
        logger.debug(f"System total RAM: {system_total_ram / giga_byte: .2f} GB")
        logger.debug(
            f"Current process RAM: {current_process_ram_usage / giga_byte: .2f} GB"
        )
        logger.debug(f"RAM reserved for OS: {ram_reserved_for_os / giga_byte: .2f} GB")
        logger.debug(f"Memory budget: {memory_budget / giga_byte: .2f} GB")
        return memory_budget

    @staticmethod
    def _get_memory_budget():
        # Check if the environment variables are set
        num_cpus = os.environ.get(TREEMS2_NUM_CPUS)
        mem_per_cpu = os.environ.get(TREEMS2_MEM_PER_CPU)
        if num_cpus and mem_per_cpu:
            logger.debug(f"Calculating memory budget based on environment variables.")
            giga_byte = 1024**3
            memory_budget = int(num_cpus) * int(mem_per_cpu) * giga_byte
            memory_budget = memory_budget / 2.0
            logger.debug(f"Memory budget: {memory_budget / giga_byte: .2f} GB")
        else:
            logger.debug(f"Estimating memory budget...")
            memory_budget = VectorStoreIndex._estimate_memory_budget()
        return memory_budget

    @staticmethod
    def _create_factory_string(vector_count: int, vector_dim: int) -> Tuple[str, int]:
        # determine the type of index (and number of clusters if applicable)
        if vector_count < 10_000:  # N < 10k
            nlist = 0
            factory_string = "IDMap"
        elif vector_count < 10**6:  # N < 1M
            nlist = min(
                math.floor(16 * math.sqrt(vector_count)), math.floor(vector_count / 39)
            )  # need a minimum of 39 training points per cluster
            factory_string = f"IVF{nlist}"
        elif vector_count < 10**7:  # N < 10M
            nlist = min(
                math.floor(20 * math.sqrt(vector_count)), math.floor(vector_count / 39)
            )
            factory_string = f"IVF{nlist}_HNSW32"
        elif vector_count < 10**8:  # N < 100M
            nlist = min(
                math.floor(26 * math.sqrt(vector_count)), math.floor(vector_count / 39)
            )
            factory_string = f"IVF{nlist}_HNSW32"
        elif vector_count <= 10**9:  # N <= 1B
            nlist = min(
                math.floor(33 * math.sqrt(vector_count)), math.floor(vector_count / 39)
            )
            factory_string = f"IVF{nlist}_HNSW32"
        else:  # N > 1B
            raise IndexConstructionError(
                f"An index can be constructed for a maximum of 1B vectors, but got {vector_count} vectors."
            )

        # determine compression
        if vector_count >= 10**6:
            memory_budget = VectorStoreIndex._get_memory_budget()
            memory_budget_per_vector = math.floor(memory_budget / vector_count)
            logger.debug(f"Memory budget per vector: {memory_budget_per_vector: .2f} B")
            m = max(memory_budget_per_vector - 16, 0)
            if m < 64:
                raise IndexConstructionError(
                    f"Not enough memory available for indexing {vector_dim} vectors. Index requires at least {(64 + 16) * vector_count} bytes, but only an estimated {memory_budget} bytes of memory available."
                )
            if m > vector_dim:
                if m > 4 * vector_dim:
                    factory_string += ",Flat"
                elif m > 2 * vector_dim:
                    factory_string += ",SQfp16"
                elif m >= vector_dim:
                    factory_string += ",SQ8"
            else:
                if 4 * m <= vector_dim:
                    factory_string = f"OPQ{m}_{4 * m}," + factory_string + f",PQ{m}"
                else:
                    factory_string = f"OPQ{m}," + factory_string + f",PQ{m}"
        else:
            factory_string += ",Flat"
        return factory_string, nlist

    def _load_training_data(self) -> Tuple[int, np.ndarray]:
        nr_of_training_points = min(39 * self.nlist, self.vector_count)

        logger.info(f"Loading {nr_of_training_points} training points from disk...")
        load_training_data_time_start = time.time()
        tracemalloc.start()

        training_data = self.vector_store.parallel_sample(nr_of_training_points)

        current, peak = tracemalloc.get_traced_memory()
        logger.debug(
            f"Loaded {nr_of_training_points} from disk, taking up {current / 1e6:.2f} MB in memory. Peak memory usage: {peak / 1e6:.2f} MB"
        )
        logger.info(
            f"Loaded all {nr_of_training_points} training points from disk in {format_execution_time(time.time() - load_training_data_time_start)}"
        )
        return nr_of_training_points, training_data

    def _train_cpu(self):
        nr_of_training_points, training_data = self._load_training_data()
        logger.info(
            f"Training index on CPU using {nr_of_training_points} training points."
        )
        train_time_start = time.time()
        self.index.train(training_data)
        logger.info(
            f"Finished training index in {format_execution_time(time.time() - train_time_start)}"
        )

    def _train_gpu(self):
        index_ivf = faiss.extract_index_ivf(self.index)
        clustering_index = faiss.index_cpu_to_all_gpus(faiss.IndexFlatIP(index_ivf.d))
        index_ivf.clustering_index = clustering_index

        nr_of_training_points, training_data = self._load_training_data()
        logger.info(
            f"Training index on GPU using {nr_of_training_points} training points."
        )
        train_time_start = time.time()
        self.index.train(training_data)
        logger.info(
            f"Finished training index in {format_execution_time(time.time() - train_time_start)}"
        )

    def _train_cpu_2level(self):
        nr_of_training_points, training_data = self._load_training_data()
        logger.info(
            f"Training index on CPU (two-level clustering) using {nr_of_training_points} training points."
        )
        train_time_start = time.time()
        clustering.train_ivf_index_with_2level(self.index, training_data)
        logger.info(
            f"Finished training index in {format_execution_time(time.time() - train_time_start)}"
        )

    def train(self, use_gpu: bool = False):
        if self.index.is_trained:
            return
        if self.nlist <= 0:
            return
        if self.vector_count < 10**7:
            self._train_cpu()
        else:
            if use_gpu:
                if faiss.get_num_gpus():
                    self._train_gpu()
                else:
                    logger.warning(
                        "No GPU found, using CPU (two-level clustering) instead for training."
                    )
                    self._train_cpu()
            else:
                self._train_cpu()

    def add(self, batch_size: int):
        """
        Add vectors present in the vector store to the FAISS index in batch.
        :param batch_size: the number of vectors in a batch
        :return:
        """
        # currently used memory
        process = psutil.Process(os.getpid())
        base_memory_usage = process.memory_info().rss

        add_time_start = time.time()

        # add vectors to the index
        with tqdm(
            desc="Vectors added to index",
            unit=f" vector",
            total=self.vector_count,
        ) as pbar:
            for vectors, ids, nr_vectors in self.vector_store.to_vector_batches(
                batch_size=batch_size
            ):
                self.index.add_with_ids(vectors, ids)
                pbar.update(nr_vectors)

        logger.info(
            f"Added all vectors to the index in {format_execution_time(time.time() - add_time_start)}"
        )

        # size of index in memory
        index_memory_usage_mb = (
            process.memory_info().rss - base_memory_usage
        ) / 1024**2
        logger.debug(f"Index takes up {index_memory_usage_mb:.2f} MB in memory.")

    def save_index(self, path):
        """
        Save the FAISS index to a file.
        :param path: filepath (str): Path to save the index.
        :return:
        """
        faiss.write_index(self.index, path)

    @staticmethod
    def load(path: str, vector_store: VectorStore) -> Optional["VectorStoreIndex"]:
        """
        Load a FAISS index from a specified file.
        """
        if not os.path.exists(path):
            return None
        try:
            vector_store_index = VectorStoreIndex(vector_store)
            index = faiss.read_index(path)  # Attempt to load the index
            vector_store_index.index = index
            return vector_store_index
        except Exception:  # Catch all exceptions silently
            return None  # Return None if loading fails

    def knn_search(
        self, k: int, nprobe: int, batch_size: int
    ) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Perform a KNN search on the FAISS index for every vector in the vector store in batches.

        :param k: (int) Number of nearest neighbors to retrieve.
        :param nprobe: (int) Number of cells (clusters) to probe (for IVF indices).
        :param batch_size: (int) Number of query vectors to process per batch.

        :yield: A tuple containing:
            - d (np.ndarray): Array of distances or similarity scores (shape: [batch_size, k]).
            - i (np.ndarray): Array of matched vector indices (shape: [batch_size, k]).
            - query_ids (np.ndarray): Array of query vector IDs (shape: [batch_size]).
        """
        self.index.nprobe = nprobe
        query_time_start = time.time()
        with tqdm(
            desc="Vectors queried", unit=" vector", total=self.vector_count
        ) as pbar:
            for (
                query_vectors,
                query_ids,
                nr_vectors,
            ) in self.vector_store.to_vector_batches(batch_size=batch_size):
                d, i = self.index.search(query_vectors, k)  # FAISS standard KNN search
                yield d, i, query_ids
                pbar.update(nr_vectors)
        logger.info(
            f"Queried all vectors in {format_execution_time(time.time() - query_time_start)}"
        )
