# Chemistry and Kit Version

> [!NOTE] 
> Only **CeleScope v2.3.0 or higher** supports analysis of **GEXSCOPE-V3** data.
> 
> Only **CeleScope v2.9.0 or higher** supports analysis of **flv_rna-V2** and **flv-V2** data.

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
