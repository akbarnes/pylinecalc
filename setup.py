import os.path
from setuptools import setup

thisdir = os.path.abspath(os.path.dirname(__file__))
version = open(os.path.join(thisdir, 'pylinecalc', 'VERSION')).read().strip()


def readme():
    with open("README.md", 'r', encoding='UTF-8') as f:
        return f.read()


setup(
    name="carsons",
    version=version,
    packages=["carsons"],
    package_data={
        '': ['VERSION'],
        'carsons': ['py.typed'],
    },
    description="A python tool for calculating overhead and underground line imedances",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    long_description=readme(),
    long_description_content_type='text/markdown',
    author="Arthur K. Barnes",
    author_email="abarnes@lanl.gov"
    url="https://github.com/bluejuniper/pylinecalc",
    keywords=["carsons", "cables", "lines", "power systems"],
    license="MIT",
    install_requires=[
        'numpy>=1.17.2',
        'PySimpleGUI>=4.16.0',
        'toml>=0.10.0',
        'jupyterlab>=2.0.1',
    ],
    zip_safe=False,
    extras_require={
        "test": [
            "pytest>=3.6",
            "pytest-cov",
            "pytest-mypy",
            "pint",
        ],
    },
)
