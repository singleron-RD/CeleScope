import setuptools
from celescope.__init__ import __VERSION__, ASSAY_LIST

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

entrys = ['celescope=celescope.celescope:main',]
for assay in ASSAY_LIST:
    entrys.append(f'multi_{assay}=celescope.{assay}.multi_{assay}:main')
entry_dict = {
        'console_scripts': entrys,
}


setuptools.setup(
    name="celescope",
    version=__VERSION__,
    author="zhouyiqi",
    author_email="zhouyiqi@singleronbio.com",
    description="Single Cell Analysis Pipelines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/singleron-RD/CeleScope",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=install_requires,
)
