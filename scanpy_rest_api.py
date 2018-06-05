import functools
import json
import multiprocessing

import numpy as np
import pandas as pd
import pandas.api.types

import scanpy.api as sc


ADATA = None
ALL_METADATA = None

def initialize(h5ad_path):
    """Set up the global ADATA object. Also parse the metadata upfront so we don't
    have to use a slow pandas function for each request.
    """

    global ADATA
    global ALL_METADATA
    ADATA = sc.read_h5ad(h5ad_path)
    ALL_METADATA = parse_metadata(ADATA)

class Float32JSONEncoder(json.JSONEncoder):
    """Python won't convert a few numpy types to JSON, so we need a different
    encoder class.
    """
    def default(self, obj):
        if isinstance(obj, np.float32):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)

def _metadata_fields(_adata):
    return _adata.obs.columns

def _is_metadata_categorical(_adata, field):
    metadata_type_dict = _adata.obs.dtypes.to_dict()
    return pandas.api.types.is_categorical_dtype(metadata_type_dict[field])

def _is_metadata_continuous(_adata, field):
    return not _is_metadata_categorical(_adata, field)

def _metadata_type(_adata, field):
    metadata_type_dict = _adata.obs.dtypes.to_dict()
    dtype = metadata_type_dict[field]
    if pandas.api.types.is_string_dtype(dtype):
        return "string"
    elif pandas.api.types.is_float_dtype(dtype):
        return "float"
    elif pandas.api.types.is_integer_dtype(dtype):
        return "int"
    return "unknown"

def parse_metadata(_adata):
    """
    Get dictionary representation of metadata for a list of cells (or all cells in
    data set if cells not set)

    :param cells: list of cells (False for all cells)
    :return: {metadata_key1: value1, metadata_key2: value1, ...}
    """

    metadata = _adata.obs.to_dict(orient="records")
    for i, idx in enumerate(_adata.obs.index):
        metadata[i]["CellName"] = idx
    return metadata

def subset_metadata(_adata):
    cell_ids = set(_adata.obs_names)
    return [m for m in ALL_METADATA if m["CellName"] in cell_ids]

def query_adata(_adata, query_string):
    """Return an adata subsetted by the query string."""
    if not query_string:
        return _adata
    # Start with all true
    cell_idx = np.ones((len(all_cells(_adata)),), dtype=bool)

    for key in query_string:
        if key in ["_nograph", "_includeisolated", "_noexpression"]:
            continue
        if key not in _metadata_fields(_adata):
            raise RuntimeError("Error: key {} not in metadata schema".format(key))

        value = query_string.getlist(key)

        if _is_metadata_categorical(_adata, key):
            key_idx = np.in1d(getattr(_adata.obs, key), value)
            cell_idx = np.logical_and(cell_idx, key_idx)
        else:
            value = value[0]
            try:
                min_, max_ = value.split(",")
            except ValueError:
                raise RuntimeError(
                    "Error: min,max format required for range for key {}, got {}".format(key, value))
            if min_ != "*":
                key_idx = np.array((getattr(_adata.obs, key) >= min_).data)
                cell_idx = np.logical_and(cell_idx, key_idx)
            if max_ != "*":
                key_idx = np.array((getattr(_adata.obs, key) <= min_).data)
                cell_idx = np.logical_and(cell_idx, key_idx)

    return _adata[cell_idx, :]

def get_metadata_ranges(_adata):
    """
    Parses through all metadata to get available values in current dataset
    :param metadata: dictionary of metadata
    :return: dictionary {categorical_key: { "options": [list of options]}, continous_key:
    {"range": {"min": min_val, "max": max_val}}}
    """
    metadata_ranges = {}

    for field in _metadata_fields(_adata):
        if _is_metadata_categorical(_adata, field):
            metadata_ranges[field] = {"options": _adata.obs.groupby(field).size().to_dict()}
        else:
            metadata_ranges[field] = {
                "range": {"min": _adata.obs[field].min(),
                          "max": _adata.obs[field].max()}}

    metadata_ranges["CellName"] = {"options": {i: 1 for i in _adata.obs.index}}

    return metadata_ranges

def create_graph(_adata, graph_method="umap", recalculate=False):

    # Run the graph method
    if recalculate:
        getattr(sc.tl, graph_method)(_adata)

    # scanpy stores the coordinates in one of the metadata tables
    df_raw = pd.DataFrame(data=_adata.obsm["X_" + graph_method], index=_adata.obs_names)

    # cellxgene wants values between 0 and 1 I think
    df_norm = (df_raw - df_raw.min()) / (df_raw.max() - df_raw.min())
    df_norm = df_norm.round(4)

    return np.hstack((df_norm.index.values.reshape(df_norm.index.shape[0], 1), df_norm.values)).tolist()

def _jsonify_list_chunk(start_end, list_to_chunk):
    start, end = start_end
    return json.dumps(list_to_chunk[start:end], cls=Float32JSONEncoder)

def jsonified_list_generator(list_):

    jsonify_list_chunk = functools.partial(_jsonify_list_chunk, list_to_chunk=list_)

    yield '['

    starts = range(0, len(list_), 20000)
    ends = range(20000, len(list_) + 20000, 20000)
    bounds = list(zip(starts, ends))

    pool = multiprocessing.Pool(8)
    counter = 0
    for result in pool.imap_unordered(jsonify_list_chunk, bounds, chunksize=20):
        trimmed_result = result[1:-1] # remove the starting and ending [ ]
        counter += 1
        if counter < len(bounds):
            yield trimmed_result + ', '
        else:
            yield trimmed_result + '] '
    pool.terminate()

def generate_cells_get(args):

    # An amazing TTFB improvement
    yield '{"data": {"reactive": true '

    nograph = bool(args.get("_nograph"))
    includeisolated = bool(args.get("_includeisolated"))

    queried_adata = query_adata(ADATA, args)

    if nograph:
        yield ', "graph": null'
    else:
        graph = create_graph(queried_adata)
        yield ', "graph": '
        for list_chunk in jsonified_list_generator(graph):
            yield list_chunk

    metadata = subset_metadata(queried_adata)

    yield ', "metadata": '
    for list_chunk in jsonified_list_generator(metadata):
        yield list_chunk

    ranges = get_metadata_ranges(queried_adata)
    yield ', "ranges": ' + json.dumps(ranges, cls=Float32JSONEncoder)

    yield ', "cellids": '
    for list_chunk in jsonified_list_generator(all_cells(queried_adata)):
        yield list_chunk

    yield ', "cellcount": ' + str(len(queried_adata.obs_names)) + '}}'

def get_schema(_adata):

    schema = {}

    for field in _metadata_fields(_adata):
        schema_dict = {
            "displayname": field,
            "variabletype": "categorical" if _is_metadata_categorical(_adata, field) else "continuous",
            "type": _metadata_type(_adata, field),
            "include": True
        }
        schema[field] = schema_dict

    schema["CellName"] = {
        "displayname": "Name",
        "variabletype": "categorical",
        "type": "string",
        "include": True
    }

    return schema

def generate_initialize_get():

    options = get_metadata_ranges(ADATA)
    genes = all_genes(ADATA)
    yield '{"data": {"reactivelimit": 999999999 '
    yield ', "ranges": ' + json.dumps(options, cls=Float32JSONEncoder)
    yield ', "cellcount": ' + str(len(all_cells(ADATA)))
    yield ', "genes": '
    for list_chunk in jsonified_list_generator(genes):
        yield list_chunk
    yield ', "schema": ' + json.dumps(get_schema(ADATA)) + '}}'


def all_genes(_adata):
    """
    Get list of all genenames in dataset
    :return: list of gene names
    """
    return _adata.var.index.tolist()


def all_cells(_adata):
    """
    Get list of all cell names in dataset
    :return: list of cell names
    """
    return _adata.obs.index.tolist()

def get_expression(_adata, cells, genes=()):
    """
    Get matrix of expression data
    :param cells: list of cell names, or all cells if empty
    :param genes: list of genes, or all genes if empty
    :return: numpy expression matrix
    """

    if cells:
        cells_idx = np.in1d(_adata.obs_names, cells)
    else:
        cells_idx = np.ones((len(all_cells(_adata)),), dtype=bool)

    if genes:
        genes_idx = np.in1d(_adata.var_names, genes)
    else:
        genes_idx = np.ones((len(all_genes(_adata)),), dtype=bool)

    subsetted_adata = _adata[cells_idx, :][:, genes_idx]

    return subsetted_adata

def remove_unexpressed_genes(expression, genes):
    """
    Filter out genes that have no expression data from expression matrix
    :param expression: numpy matrix will expression data
    :param genes: list of genes names
    :return: tuple, filtered expression numpy matrix, list of genes
    """
    genes_expressed = expression.any(axis=0)
    genes = [genes[idx] for idx, val in enumerate(genes_expressed) if val]
    return expression[:, genes_expressed], genes

def parse_exp_data(_adata, cells=(), genes=(), limit=0, unexpressed_genes=False):
    """
    Get expression data for set of cells and genes
    :param cells: list of cell names
    :param genes: list of gene names
    :param limit: optional, limit number of genes returned
    :param unexpressed_genes: boolean, filter out genes with no expression
    :return: json: genes: list of genes, cells: expression matrix:,
        nonzero_gene_count: number of expressed genes
    """

    expression = get_expression(_adata, cells, genes).X

    if len(expression.shape) == 1:
        expression = expression.reshape(expression.shape[0], 1)

    if not genes:
        genes = all_genes(_adata)
    if not cells:
        cells = all_cells(_adata)

    if not unexpressed_genes:
        expression, genes = remove_unexpressed_genes(expression, genes)
    if limit and len(genes) > limit:
        genes = genes[:limit]
        expression = expression[:, :limit]
    cell_data = []
    for idx, cell in enumerate(cells):
        cell_data.append({
            "cellname": cell,
            "e": list(expression[idx]),
        })
    return {
        "genes": genes,
        "cells": cell_data,
        "nonzero_gene_count": int(np.sum(expression.any(axis=1)))
    }
