from setuptools import setup

setup (
    name = 'pydiskvert',
    version = '170611',
    author = 'Dominik Gronkiewicz',
    author_email = 'gronki@gmail.com',
    description = u"Calculate vertical structure of accretion disks",
    license = "MIT",
    packages = [ 'diskvert' ],
    scripts = [
        'scripts/col2python',
        'scripts/diskvert-cooling2D',
        'scripts/diskvert-new-plot',
    ],
    install_requires = [
        'numpy',
    ],
)
