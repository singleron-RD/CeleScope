import setuptools
from celescope.__init__ import __VERSION__, ASSAY_DICT
with open("README.md", "r") as fh:
    long_description = fh.read()

entrys = ['celescope=celescope.celescope:main',]
for assay in ASSAY_DICT:
    entrys.append(f'multi_{assay}=celescope.{assay}.multi_{assay}:main')
entry_dict = {
        'console_scripts': entrys,
}


setuptools.setup(
    name="celescope",
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
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=[
        'cutadapt==1.17',
        'pysam>=0.15.0',
        'scipy>=0.19.1',
        'numpy>=1.15.4',
        'pandas>=0.23.4',
        'jinja2>=2.10',
        'matplotlib>=2.2.2',
        'xopen>=0.5.0',
        'editdistance>=0.5.3'
    ]
)
