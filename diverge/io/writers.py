"""
Functions for writing analysis results to various file formats.
"""
from typing import Optional, Dict, Any, Union
import os
import pandas as pd
import json

def save_results_to_csv(
    df: pd.DataFrame,
    filename: str,
    index: bool = True,
    **kwargs
) -> str:
    """
    Save pandas DataFrame to CSV file.
    
    Args:
        df: DataFrame to save
        filename: Output filename
        index: Whether to write row names (index)
        **kwargs: Additional arguments passed to pandas.DataFrame.to_csv()
        
    Returns:
        Path to saved file
    """
    # 确保文件名有正确的扩展名
    if not filename.endswith('.csv'):
        filename += '.csv'
        
    # 确保目录存在
    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    
    # 保存数据
    df.to_csv(filename, index=index, **kwargs)
    return filename

def save_results_to_excel(
    result_dict: Dict[str, pd.DataFrame],
    filename: str,
    **kwargs
) -> str:
    """
    Save multiple DataFrames to Excel file with multiple sheets.
    
    Args:
        result_dict: Dictionary mapping sheet names to DataFrames
        filename: Output filename
        **kwargs: Additional arguments passed to pandas.ExcelWriter
        
    Returns:
        Path to saved file
    """
    # 确保文件名有正确的扩展名
    if not filename.endswith(('.xlsx', '.xls')):
        filename += '.xlsx'
        
    # 确保目录存在
    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    
    # 保存数据到多个工作表
    with pd.ExcelWriter(filename, **kwargs) as writer:
        for sheet_name, df in result_dict.items():
            df.to_excel(writer, sheet_name=sheet_name)
    
    return filename 