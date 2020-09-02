from setuptools import setup, find_packages


setup(

    name='varcomb',
    version='0.0.1',

    packages=find_packages('src'),
    package_dir={'': 'src'},

    test_suite='tests',

    entry_points={
        'console_scripts': ['somaseq = somaseq.client:run']
    },

    python_requires='>=3.8',

    install_requires=[
        'click'
    ],

    author='Daniel Thaagaard Andreasen',
    author_email='daniel@clin.au.dk',
    license='MIT'

)
