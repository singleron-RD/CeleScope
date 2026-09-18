# Chemistry and Kit Version

> [!NOTE] 
> Analysis of scRNA **GEXSCOPE-V3** is supported only in **CeleScope v2.3.0 or later**.
>
> Analysis of scRNA combined with full-length VDJ **flv_rna-V2** and **flv-V2** is supported only in **CeleScope v2.9.0 or later**.


Upgrade command:  
```bash
mamba activate celescope # or source activate celescope
pip install --upgrade celescope
```

## Mapping of Kit Version to Chemistry  

Starting from CeleScope v2.3.0, the **chemistry** field in the HTML report is aligned with the **Kit Version**.  

| Kit Version | CeleScope >= v2.3.0  | CeleScope < v2.3.0 |
|------------|----------------------|--------------------|
| **Microbead** | `GEXSCOPE-MicroBead` | `scopeV1`          |
| **V1**       | `GEXSCOPE-V1`        | `scopeV2*`         |
| **V2**       | `GEXSCOPE-V2`        | `scopeV3.0.1`      |
| **V3**       | `GEXSCOPE-V3`        | *Not supported*    |

## Chemistry Patterns  

Chemistry patterns can be found [here](https://github.com/singleron-RD/CeleScope/blob/master/celescope/chemistry_dict.py).  

- **`C`**: Cell Barcode  
- **`U`**: Unique Molecular Identifier (UMI)  
- **`L`**: Linker (a fixed sequence separating multiple barcode segments)  

For example, in the pattern:  
```plaintext
"GEXSCOPE-V1": "C8L16C8L16C8L1U12"
```
- The first **8 bp** of R1 is the **cell barcode**.  
- The next **16 bp** is a **linker**.  
- Then, another **8 bp** of **cell barcode**, and so on.  

> [!NOTE] 
> To maintain sequencing base balance, `GEXSCOPE-V3` and `flv_rna-V2` includes an additional **0–3 bp** sequence before the first barcode segment.

## GEXSCOPE-V3 Structure

### Sequence Pattern

```plaintext
C9L6C9L6C9L1U12
```

Due to the **0–3 bp staggered offset** at the beginning of Read 1 (for base balance and sequence diversity), the actual pattern in Read 1 is:

```plaintext
[0-3 bp offset] + C9 + L6 + C9 + L6 + C9 + L1(1bp C) + U12
```

### Region Breakdown

| Region | Pattern | Length | Description | Whitelist |
|--------|---------|--------|-------------|-----------|
| Initial Offset | — | 0–3 bp | Staggered random nucleotides at the start of R1 | — |
| Barcode 1 (C1) | C9 | 9 bp | First cell barcode segment | [bc1.txt](../celescope/data/chemistry/GEXSCOPE-V3/bc1.txt) |
| Linker 1 (L1) | L6 | 6 bp | Fixed sequence: `ACGATG` | [linker1.txt](../celescope/data/chemistry/GEXSCOPE-V3/linker1.txt) |
| Barcode 2 (C2) | C9 | 9 bp | Second cell barcode segment | [bc2.txt](../celescope/data/chemistry/GEXSCOPE-V3/bc2.txt) |
| Linker 2 (L2) | L6 | 6 bp | Fixed sequence: `CATAGT` | [linker2.txt](../celescope/data/chemistry/GEXSCOPE-V3/linker2.txt) |
| Barcode 3 (C3) | C9 | 9 bp | Third cell barcode segment | [bc3.txt](../celescope/data/chemistry/GEXSCOPE-V3/bc3.txt) |
| Spacer | L1(1bp C) | 1 bp | Spacer nucleotide following Barcode 3 | — |
| UMI | U12 | 12 bp | Unique Molecular Identifier for transcript deduplication | — |

### Key Specifications

- **Total Cell Barcode**: 27 bp across three 9 bp segments (`C9 + C9 + C9`)
- **UMI**: 12 bp fixed length
- **Dynamic Location Handling**: Because of the 0–3 bp offset, tools like STARsolo process GEXSCOPE-V3 under `CB_UMI_Complex` mode using `soloAdapterSequence` (`NNNNNNNNNACGATGNNNNNNNNNCATAGT`) to anchor and locate the barcode segments relative to the linkers