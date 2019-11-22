from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name = "PESViewer",
    version = "1.0",

    packages = find_packages(),
    entry_points={'console_scripts':[
        'pesviewer = pesviewer.pesviewer:main',
        ]},
    install_requires=['matplotlib','numpy','Pillow'],

    author="Ruben Van de Vijver",
    author_email = "vandevijver.ruben@gmail.com",

    long_description=long_description,
    long_description_type='text/markdown',
    description = "Depiction of Potential Energy Surfaces.",
    license = "MIT",
    url = "https://github.com/rubenvdvijver/PESViewer.git",
)
