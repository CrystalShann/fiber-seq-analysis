import sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Defaults preserve the original behavior (the 100 bp accessibility file).
# CLI usage: convert_accessibility_to_parquet.py <src_tsv_gz> <dst_parquet>
DEFAULT_SRC = "/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_quartiles/AL10_bc2178_19130_1kb_bin10_modthresh0.9_all/tables/AL10_bc2178_19130_tss_read_accessibility.tsv.gz"
DEFAULT_DST = "/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_quartiles/AL10_bc2178_19130_1kb_bin10_modthresh0.9_all/tables/AL10_bc2178_19130_tss_read_accessibility.parquet"

if len(sys.argv) == 1:
    src, dst = DEFAULT_SRC, DEFAULT_DST
elif len(sys.argv) == 3:
    src, dst = sys.argv[1], sys.argv[2]
else:
    sys.exit("usage: convert_accessibility_to_parquet.py [<src_tsv_gz> <dst_parquet>]")

dtypes = {
    "read_id": "string",
    "tss_uid": "category",
    "gene_id": "category",
    "gene_id_base": "category",
    "gene_name": "category",
    "chromosome": "category",
    "tss_coordinate": "int32",
    "tss_strand": "category",
    "tss_source": "category",
    "PI": "float32",
    "pausing_group": "category",
    "n_positions_covered": "int32",
    "n_m6a_modified": "int32",
    "n_m6a_unmodified": "int32",
    "tss_accessible": "boolean",
}

writer = None
total_rows = 0
for i, chunk in enumerate(pd.read_csv(src, sep="\t", dtype=dtypes, chunksize=2_000_000)):
    table = pa.Table.from_pandas(chunk, preserve_index=False)
    if writer is None:
        writer = pq.ParquetWriter(dst, table.schema, compression="zstd")
    writer.write_table(table)
    total_rows += len(chunk)
    print(f"chunk {i}: {len(chunk):,} rows written (total {total_rows:,})", flush=True)

writer.close()
print(f"done: wrote {total_rows:,} rows from {src} to {dst}", flush=True)
