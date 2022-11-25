from setuptools import setup, find_packages

setup(
    name="ARstack",
    version="0.0.1",
    author="ARStack",
    description="Stacking model for protein structure prediction",
    url="https://github.com/HangL-39/ARStack",
    packages=find_packages(),
    python_requires='>=3.7',
    install_requires=['Bio == 1.5.1',
                      'tensorflow == 2.11.0',
                      'keras == 2.11.0',
                      'numpy == 1.23.5',
                      'pandas == 1.5.2',
                      'scikit_learn == 1.1.3',
                      'scipy == 1.9.3',
                      'seaborn == 0.12.1',
                      'setuptools == 65.5.0',
                      ],
)
