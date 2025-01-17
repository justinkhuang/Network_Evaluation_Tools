"""
Setup module adapted from setuptools code. See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

setup(
	name='network_evaluation_tools',
	version='1.0.2',
	description='Module to perform patient and molecular network evaluation as described in Huang and Carlin, et al. 2018',
	url='https://github.com/idekerlab/Network_Evaluation_Tools',
	author='Justin Huang',
	author_email='jkh013@ucsd.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Software Development :: Build Tools',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2.7'
	],
	packages=find_packages(exclude=['copy', 'itertools', 'os', 're', 'time']),
	install_requires=[
        'argparse>=1.1',
        'networkx>=2.1',
        'numpy>=1.16.3',
        'matplotlib>=3.0.3',
        'pandas>=0.24.2',
        'requests>=2.21.0',
        'scipy>=1.2.1',
        'scikit-learn>=0.20.3',
        'seaborn>=0.7.1']
)