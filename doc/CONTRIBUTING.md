## Pull Requests
Create pull requests to `dev` branch

## Lint code
Before pull requests, you should lint your code with the following command:
```
pip install pylint
# W1618 (no-absolute-import)
# E1101 (no-member)
# W1633 (round-builtin)
# W1619 (old-division)
# W0105 (String statement has no effect)
# W0511 TODO!
# E1130 bad operand type for unary ~: _isnan (invalid-unary-operand-type)
# W0212 protected-access
# W0221 arguments-differ
pylint --disable=all --enable=E,W --disable=W1618,E1101,W1633,W1619,W0105,W0511,E1130,W0212,W0221 --jobs=2 celescope
```
Your code should be rated at 10(i.e. no error or warning). 

## Add a new assay

1. Add the assay in `celescope.__init__.ASSAY_LIST` . If the assay is ready to be released, add it to `RELEASED_ASSAYS`.
2. Create a submodule under celescope with the exact assay name.
3. Add new steps under this submodule(folder).

## Add a new step for an assay
When you add a new step, you need to

  - Create a class which inherit from `celescope.tools.step.Step`. 

  - Create a function with the same name as the module(file name). The main function `celescope` uses this function to run each step.

  - Create a parser function with the name `get_opts_{module_name}`. `celescope` command line interface uses this function. The `sub_program` argument in this function hides all the arguments that you do not want to show in the `multi_{assay}` interface.

  - If a HTML report is needed, you need to add `{assay}/base.html` in `celescope/templates`. Take a look at `celescope/rna/base.html`. 

For example, in `celescope.tools.cutadapt`:

```
from celescope.tools.step import Step, s_common
from celescope.tools import utils


class Cutadapt(Step):
    """
    Features
    - Trim adapters in R2 reads with cutadapt. Default adapters includes:
	- polyT=A{18}, 18 A bases. 
	- p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.

    Output
    - `cutadapt.log` Cutadapt output log file.
    - `{sample}_clean_2.fq.gz` R2 reads file without adapters.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        {some init code}

    @utils.add_log
    def run(self):
        {some code to run}


@utils.add_log
def cutadapt(args):
    with Cutadapt(args, display_title="Trimming") as runner:
        runner.run()


def get_opts_cutadapt(parser, sub_program):
    parser.add_argument('--adapter_fasta', help='Addtional adapter fasta file.')
    parser.add_argument(
        '--minimum_length',
        help='Default `20`. Discard processed reads that are shorter than LENGTH.', 
        default=20
    )
    {other arguments}
    if sub_program:
        parser.add_argument('--fq', help='Required. R2 reads from step Barcode.', required=True)
        parser.add_argument('--gzip', help="Output gzipped fastq", action='store_true')
        parser = s_common(parser)
    return parser
```

## Reusable steps
All the steps under `celescope/tools` are reusable. For example, all the assays uses the `sample`, `barcode` and `cutadapt` . You don't need to copy the code into each assay. Just put them in  `celescope.{assay}.__init__.STEPS`

celescope/rna/__init__.py
```
STEPS = [
    'mkref',
    'sample',
    'barcode',
    'cutadapt',
    'star',
    "featureCounts",
    "count",
    'analysis']
__ASSAY__ = 'rna'
```

## Docs
Under CeleScope root direcotry, run

`python scripts/generate_docs.py`

This will generate documents for each step. The generated docs are in the `docs` folder. It will collect:

- Docstring of the step class. The Docstring should have sections named `Features` and `Output`.
- Help infomation in `get_opts_{module_name}`
  
Released assays will be added to `manual.md`. 


## Tests
If you add new steps, you need to create a small data for integration tests. 
 - To run sample tests, See https://github.com/singleron-RD/celescope_test_script

Then you need to create your own test based on this example. 
 - You can use `scripts/extract_read.py` to extract reads assigned to some cells. This script may help to create test data.