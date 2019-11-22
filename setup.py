from setuptools import setup, find_packages

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
    
    description = "Depiction of Potential Energy Surfaces.",
    license = "MIT",
    url = "https://github.com/rubenvdvijver/PESViewer.git",
)
