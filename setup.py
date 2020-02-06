from setuptools import setup, find_packages

setup(
    name="dmdpy",
    version="1.0.0",
    author=["Matthew R. Hennefarth", "David Reilley"],
    packages=find_packages(),
    package_data={'': ['.json', '*.j2', '*.config']},
    entry_points={
        'console_scripts' : [
            'relabelpdb.py=dmdpy.bin.relabelpdb:main',
            'setupdmd.py=dmdpy.bin.setupdmd:main',
            'rundmd.py=dmdpy.bin.rundmd:main',
            'm2p=dmdpy.bin.movietopdb:main'
        ]
    },
)
