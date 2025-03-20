"""
Parallel processing utilities for DIVERGE.
"""
from typing import List, Callable, Any, TypeVar, Generic
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import numpy as np

T = TypeVar('T')
R = TypeVar('R')

def parallel_map(
    func: Callable[[T], R],
    items: List[T],
    n_jobs: int = None,
    use_threads: bool = False
) -> List[R]:
    """
    Execute a function in parallel across multiple items.
    
    Args:
        func: Function to apply to each item
        items: List of items to process
        n_jobs: Number of parallel jobs (default: number of cores)
        use_threads: Whether to use threads instead of processes
        
    Returns:
        List of results, one per input item
    """
    if n_jobs is None:
        n_jobs = mp.cpu_count()
    
    # 选择执行器
    executor_class = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    
    with executor_class(max_workers=n_jobs) as executor:
        results = list(executor.map(func, items))
    
    return results

def chunked_parallel(
    func: Callable[[List[T]], List[R]],
    items: List[T],
    chunk_size: int = 100,
    n_jobs: int = None
) -> List[R]:
    """
    Process items in parallel chunks for better performance with large datasets.
    
    Args:
        func: Function that processes a batch of items
        items: List of all items
        chunk_size: Size of each batch
        n_jobs: Number of parallel jobs
        
    Returns:
        Combined results from all batches
    """
    # 将数据分成块
    chunks = [items[i:i+chunk_size] for i in range(0, len(items), chunk_size)]
    
    # 并行处理每个块
    chunk_results = parallel_map(func, chunks, n_jobs=n_jobs)
    
    # 合并结果
    all_results = []
    for result in chunk_results:
        all_results.extend(result)
    
    return all_results 