import setuptools
from celescope.tools.__version__ import __VERSION__
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="celescope", # Replace with your own username
    version=__VERSION__,
    author="zhouyiqi",
    author_email="zhouyiqi@singleronbio.com",
    description="GEXSCOPE Single cell analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zhouyiqi91/CeleScope",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    package_data={
        'celescope.tools': ['*.R'],     # All R files 
        '': ['templates/*'],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'celescope=celescope.celescope:main',
        ],
    },
    install_requires=[
        'cutadapt==1.17',
        'pysam==0.15.3',
        'scipy==0.19.1',
        'numpy==1.15.4',
        'pandas==0.23.4',
        'jinja2==2.10',
        'matplotlib==2.2.2',
        'xopen==0.5.0',
    ]
)
