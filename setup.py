from setuptools import find_packages, setup

with open('urbansprawl/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break
    else:
        version = '0.0.1'

with open('README.md', 'rb') as f:
    readme = f.read().decode('utf-8')

install_requires = [
    'psutil',
    'numpy<=1.14.1',
    'pandas',
    'matplotlib',
    'shapely',
    'geopandas',
    'scikit-learn',
    'tensorflow<=1.10.0',
    'keras',
    'networkx',
    'osmnx',
    'jupyter'
]

setup(
    name='urbansprawl',
    keywords=['urbansprawl', 'land use mix', 'gis', 'spatial analysis', 'machine learning', 'openstreetmap', 'population density', 'population downscaling', 'neural networks'],
    version=version,
    description='The urbansprawl project provides an open source framework for assessing urban sprawl using open data',
    long_description=readme,
    author='Luciano Gervasoni',
    author_email='gervasoni.luc@gmail.com',
    maintainer='Luciano Gervasoni',
    maintainer_email='gervasoni.luc@gmail.com',
    license='MIT',
    url='https://github.com/lgervasoni/urbansprawl',
    entry_points={ },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    install_requires=install_requires,
    # pip install -e .[dev]
    extras_require={'dev': ['pytest', 'flake8', 'ipython', 'ipdb']},
    packages=find_packages(exclude=['examples']),
)
