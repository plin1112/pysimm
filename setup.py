from distutils.core import setup

setup(
    name = 'pysimm',
    packages = [
        'pysimm',
        'pysimm.apps',
        'pysimm.apps.solvate',
        'pysimm.forcefield',
        'pysimm.models.monomers',
        'pysimm.models.monomers.gaff',
        'pysimm.models.monomers.gaff2',
        'pysimm.models.monomers.dreiding'
    ],
    version = '0.2dev',
    description = 'python simulation interface for molecular modeling',
    author = 'Michael E. Fortunato',
    author_email = 'mef231@gmail.com',
    url = 'https://github.com/polysimtools/pysimm',
    download_url = 'https://github.com/polysimtools/pysimm/archive/0.1.tar.gz',
    scripts = ['bin/pysimm'],
    package_data = {'pysimm': [
        'forcefield/dat/*.json', 
        'forcefield/dat/*.xml',
        'apps/solvate/dat/*.json'
    ]}
)
