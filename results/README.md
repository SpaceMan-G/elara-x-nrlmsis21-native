# Paper-linked model results

This directory is reserved for publication-safe model-specific outputs.

Use the controlled Elara X atmospheric repository result router rather than
manually copying complete private execution directories into this repository.

Expected hierarchy:

```text
results/
└── paper/
    └── <paper-id>/
        └── <execution-id>/
            ├── RUN_MANIFEST.json
            ├── CHECKSUMS.sha256
            ├── results.csv
            ├── summary.json
            └── figures/
```

Not every execution requires every example file. The manifest and checksums
are mandatory whenever a result set is published.
