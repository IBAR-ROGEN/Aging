# POD5 `run_infos` vs Reader API (0.3.49)

Read-only probe. Version pins were not changed.

## Provenance

| Field | Value |
|-------|-------|
| timestamp_utc | 2026-11-24T09:20:15Z |
| workspace HEAD | `430f0558a0ad899c0f7dab0ff63fc7128ee0bf71` |
| origin/main | `f7c4c3d4f59b57cfbe648012735768ebadf678cb` |
| Python | 3.12.8 (`uv run`) |
| pod5 | **0.3.49** (`importlib.metadata.version("pod5")`) |
| pod5 path | `.venv/lib/python3.12/site-packages/pod5/` |
| intake path | `src/rogen_aging/intake/leakcheck.py` `_scan_pod5` |

## Intake expectation

`_scan_pod5` does:

```python
run_infos = getattr(reader, "run_infos", None)
if run_infos is None:
    return POD5_NOT_SCANNED  # "pod5 metadata not scanned"
```

## Minimal reproduce

```python
import importlib
import importlib.metadata

version = importlib.metadata.version("pod5")
module = importlib.import_module("pod5")
reader_cls = getattr(module, "Reader")

assert version == "0.3.49"
assert getattr(reader_cls, "run_infos", None) is None
assert not hasattr(reader_cls, "run_infos")

try:
    _ = reader_cls.run_infos
except AttributeError as exc:
    # Observed: type object 'Reader' has no attribute 'run_infos'
    assert "run_infos" in str(exc)
else:
    raise AssertionError("expected AttributeError for Reader.run_infos")

print("pod5_version", version)
print("intake_would_return", "pod5 metadata not scanned")
```

Observed output:

```
pod5_version=0.3.49
Reader.run_infos=None
hasattr(Reader,'run_infos')=False
direct_access_AttributeError=type object 'Reader' has no attribute 'run_infos'
intake_behavior: getattr(reader,'run_infos',None) is None -> return 'pod5 metadata not scanned'
```

No POD5 file is required to show the missing attribute; `getattr(Reader, "run_infos", None)` is already `None`.

## Available API (pod5 0.3.49)

### `pod5.Reader` public attributes

`batch_count`, `close`, `file_identifier`, `file_version`, `file_version_pre_migration`, `get_batch`, `inner_file_reader`, `is_vbz_compressed`, `logical_file_version`, `logical_read_table_version`, `num_reads`, `original_file_version`, `path`, `physical_read_table_version`, `read_batches`, `read_ids`, `read_ids_raw`, `read_table`, `reads`, `reads_table_version`, `run_info_table`, `signal_batch_row_count`, `signal_table`, `writing_software`

- Present: `run_info_table`, `reads`
- Absent: `run_infos`, `run_info`, `get_run_info`

### `pod5.ReadRecord`

- Present: `run_info`, `run_info_index` (per-read run metadata; used by tests that iterate `reader.reads()`)

## Conclusion

On the pinned pod5 **0.3.49**, `Reader` has no `run_infos`. Intake therefore returns `pod5 metadata not scanned` whenever POD5 is installed, without scanning run metadata via that attribute. Per-read `ReadRecord.run_info` and `Reader.run_info_table` exist as alternate surfaces; intake does not use them today.
