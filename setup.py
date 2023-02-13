from setuptools import setup, find_packages
from import_serial import __version__

with open('README.rst') as f:
    long_description = ''.join(f.readlines())


setup(
    name='import_serial',
    version=__version__,
    description='Import diffractio data from serial macromolecular crystallography',
    long_description=long_description,
    author='Martin Maly',
    author_email='martin.maly@soton.ac.uk',
    keywords='macromolecular crystallography, research',
    license='GNU Lesser General Public License v3 (LGPLv3)',
    url='https://www.ccp4.ac.uk/',
    packages=find_packages(),
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'import_serial = import_serial.import_serial:run',
        ]
    },
    # install_requires=['numpy', 'matplotlib'],
    # package_data={'pairef': ['static/*.css']},
    # data_files=[('bitmaps', ['pairef/static/pairef_logo_64.png'])],
    # include_package_data=True,
)
