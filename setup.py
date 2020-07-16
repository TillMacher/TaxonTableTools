import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="taxontabletools", # Replace with your own username
    version="1.0.5",
    author="Till-Hendrik Macher",
    author_email="till-hendrik.macher@uni-due.de",
    description="taxontabletools: A comprehensive and user-friendly program to analyse and visualise DNA metabarcoding data",
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
                        'openpyxl>=3.0.3',
                        'xlsxwriter>=1.2.7',
                        'biopython',
                        'scikit-bio'],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
