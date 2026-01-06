## What if CeleScope reports an error?

When CeleScope reports an error, using GitHub Issues can be an effective way to solve or report the error.

Before creating a new issue, search the existing issues with keywords(e.g., error message) to see if someone else has already reported the same problem. 

## How to install a previous version instead of the latest version?

If Docker is available on your machine, using the Docker image is the most convenient option:

https://quay.io/repository/singleron-rd/celescope?tab=tags

An example of running celescope interactively with Docker is shown below:
```bash 
docker run -it --rm \
  --user $(id -u):$(id -g) \
  -v /SGRNJ06/randd/:/SGRNJ06/randd/ \
  quay.io/singleron-rd/celescope:v2.10.3 /bin/bash
```

Alternatively, mamba/conda can be used.

1. Download an earlier version of `conda_pkgs.txt`: just change `master` to the version number(e.g. v1.15.4) in the link to that version 
```
wget https://raw.githubusercontent.com/singleron-RD/CeleScope/v1.15.4/conda_pkgs.txt
# create the env
mamba create -n celescope1.15.4 -y --file conda_pkgs.txt
```

2. Install celescope
```
mamba activate celescope1.15.4
pip install celescope==1.15.4
```

## "auto chemistry detection failed"

This issue is usually caused by using an outdated version of celescope, which is unable to detect the latest [GEXSCOPE-V3 chemistry](./chemistry.md).

If you are currently using celescope v1, you can install version 1.15.4 in the existing conda environment:
`pip install celescope==1.15.4`

However, it is strongly recommended to install the latest version of celescope.
celescope v2 uses STARsolo, which is more than 5× faster than v1 and requires significantly less disk space and memory.

Because celescope v2 upgrades the STAR version, it requires a new conda environment and genome build, and therefore cannot be installed via pip in a celescope v1 conda environment.

## "Failed building wheel for {package_name}" when `pip install celescope`

This error occurs because some python package contains C-extensions that your system is trying to compile from source, but it lacks the necessary C++ compiler or development headers. To resolve this without setting up a complex build environment, you should install the pre-compiled binary package by running:
```Bash
pip install {package_name} --only-binary :all:
```

Using the --only-binary :all: flag forces pip to download a version that is already "built," bypassing the compilation step entirely and ensuring a smooth installation.

## Create .loom file from celescope output

You can use the Python code in the gist below to create a .loom file for analysis with [Veloctyo](https://velocyto.org/):

https://gist.github.com/zhouyiqi91/8535eed2c8fa8b415d252163efde7a18

## Doublet rate

For standard Singleron chips (cell number ≤ 20k), the doublet rate can be estimated as:

`captured cell number / 1k × 0.4%`

For example, if 10k cells are captured, the estimated doublet rate is:

`10k / 1k × 0.4% = 4%`
