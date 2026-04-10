# multiqc-celescope

When you have many samples and want a unified summary of analysis metrics, install `multiqc-celescope` and generate a combined HTML report.

## Install

```bash
pip install multiqc-celescope --only-binary=:all:
```

## Generate HTML report

Run `multiqc` with the `celescope` module enabled and point it to the analysis directory:

```bash
multiqc -m celescope -v {path}
```

Replace `{path}` with the root directory containing your sample output folders. This command scans all `celescope` analysis `.data.json` files under the given path and collects sample metrics from them. It will create a `multiqc_report.html` file summarizing metrics across all samples.

### Multiple analysis paths

If you have more than one analysis directory, pass multiple paths to `multiqc`:

```bash
multiqc -m celescope -v /path/to/dir1 /path/to/dir2
```

### MultiQC
For more details on MultiQC usage, please refer to:
https://github.com/MultiQC/MultiQC


