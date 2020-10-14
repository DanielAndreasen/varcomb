from setuptools import find_packages, setup

setup(

    name='varcomb',
    version='0.0.3',

    packages=find_packages('src'),
    package_dir={'': 'src'},

    test_suite='tests',

    entry_points={
        'console_scripts': ['varcomb = varcomb.client:run']
    },

    python_requires='>=3.8',

    install_requires=[
        'click'
    ],

    author='Daniel Thaagaard Andreasen',
    author_email='daniel@clin.au.dk',
    license='MIT'

)
