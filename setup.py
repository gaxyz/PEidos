import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()
    
setuptools.setup(
    name = "peidos-grijo",
    version = "0.0.1",
    author = "Gaston Rijo",
    author_email = "grijo@pasteur.edu.uy",
    description = "Python writes Eidos scripts.",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/gaxyz/PEidos",
    packages = setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"],
    python_requires = '>=3.6',    
    )