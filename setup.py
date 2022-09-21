from setuptools import setup

setup(name='myharmonizer',
      packages=['myharmonizer'],
      include_package_data=True,
      version='0.31',
      description='Myharmonizer',
      url='https://github.com/bicbioeng/myharmonizer',
      author='TD',
      author_email='tuyen.do@usd.edu',
      keywords =  'FSAS project utilities',
      license='MIT',
      python_requires='>=3.10',
      install_requires=['scikit-learn','pandas','numpy','scipy','jupyter','matplotlib','Pillow','itertools', 'seaborn','pathlib','rpy2'],
      zip_safe=False,
      classifiers=['Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3.10'])