from setuptools import setup

setup(name='xtal',
      version='0.9', # Update __init__.py if the version changes!
      description='Tools to manipulate atomic trajectories',
      author='Aravind Krishnamoorthy',
      author_email='arvk@users.noreply.github.com',
      license="MIT",
      url='https://github.com/USCCACS/xtal',
      packages=['xtal'],
      install_requires=['progressbar2'],
     )
