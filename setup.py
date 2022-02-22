import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="taxontabletools", # Replace with your own username
    version="1.4.0",
    author="Till-Hendrik Macher",
    author_email="till-hendrik.macher@uni-due.de",
    description="TaxonTableTools - A comprehensive, platform-independent graphical user interface software to explore and visualise DNA metabarcoding data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/taxontabletools",
    packages=setuptools.find_packages(),
    license = 'MIT',
    install_requires = ['pySimpleGUI>=4.15.2',
                        'pandas>=0.25.3',
                        'numpy>=1.18.1',
                        'matplotlib-venn>=0.11.5',
                        'xlrd>=1.2.0',
                        'tqdm>=4.60.0',
                        'openpyxl>=3.0.3',
                        'xlsxwriter>=1.2.7',
                        'biopython>=1.77',
                        'scikit-bio>=0.5.6',
                        'requests_html>=0.10.0',
                        'scipy>=1.5.1',
                        'shapely>=1.7.1',
                        'plotly>=4.9.0',
                        'kaleido>=0.0.3',
                        'statsmodels>=0.12.1',
                        'lastversion>=1.3.3',
                        'networkx>=2.4'],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points = {
        "console_scripts" : [
            "taxontabletools = taxontabletools.__main__:main",
        ]
    },
)
