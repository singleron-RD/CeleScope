## Pull Requests
Create pull requests to `dev` branch

## Lint code
Before pull requests, you should lint your code with the following command:
```
pip install pylint
# lint
# W1618 (no-absolute-import)
# E1101 (no-member)
# W1633 (round-builtin)
# W1619 (old-division)
# W0105 (String statement has no effect)
# W0511 TODO!
# E1130 bad operand type for unary ~: _isnan (invalid-unary-operand-type)
# W0212 Access to a protected member _option_string_actions of a client class (protected-access)
pylint --disable=all --enable=E,W --disable=W1618,E1101,W1633,W1619,W0105,W0511,E1130,W0212 --jobs=8 celescope
```
Your code should be rated at 10(i.e. no error or warning). 

## Write a new step
When you add a new step, you need to
  - Create a step class which inherit from `celescope.tools.step.Step`. 
  - Create a function with the same name of the module. The main function `celescope` uses this function to run each step.
  - Create a parser function with the name `get_opts_{module_name}`. `celescope` command line interface uses this function. The `sub_program` argument in this function hides all the arguments that you do not want to show in the `multi_{assay}` interface.

For example, in `celescope.tools.cutadapt`:
```
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

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        {some init code}

    def run(self):
        {some code to run}


def cutadapt(args):
    step_name = "cutadapt"
    cutadapt_obj = Cutadapt(args, step_name)
    cutadapt_obj.run()

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

## Docs
There is a python script at the root of this repo `generate_docs.py` to generate documents for each released step. The generated docs are in the `docs` folder. It will collect:
- Docstring of the step class. The Docstring should have sections named `Features` and `Output`.
- Help infomation in `get_opts_{module_name}`
  
Released assays will be added to `manual.md`.

## Tests
If you add new steps, you need to create a small data for integration tests. There is a test example in `celescope/tests/test_multi.py`. To run this example:

1. Copy `/SGRNJ03/randd/user/zhouyiqi/multi_tests/test_folder` to {test_dir}
2. Run `pytest`
```
Install pytest
>>> pip install pytest
Run all
>>> pytest -s celescope/tests/test_multi.py --test_dir {test_dir}
Run some tests
>>> pytest -s celescope/tests/test_multi.py --test_dir {test_dir} --assays rna,tag
```

Then you need to create your own test based on this example.