from setuptools import setup

setup (
    name = 'pydiskvert',
    version = '170623',
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
    entry_points = {
        'console_scripts': [
            'dvpl-alpha-rx=diskvert.plot_alpharx:main_plotmulti',
            'dvpl-alpha-rx-combine=diskvert.plot_alpharx:main_plotcumul',
        ],
    },
    install_requires = [
        'numpy', 'matplotlib', 'sympy', 'ipython',
    ],
)
