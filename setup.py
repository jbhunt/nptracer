from setuptools import setup, find_packages

setup(
    name='nptracer',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pynrrd',
        'numpy', 
        'allensdk'
    ],
    author='Josh Hunt',
    author_email='jbhunt303@gmail.com',
    description='A Python package for localizing units from extracellular recordings',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yjbhunt/nptracer',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
