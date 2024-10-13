from setuptools import setup, find_packages

setup(
    name='monte_carlo_simulation',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib'
    ],
    entry_points={
        'console_scripts': [
            'monte_carlo_simulation=monte_carlo_simulation.__main__:main',
        ],
    },
)
