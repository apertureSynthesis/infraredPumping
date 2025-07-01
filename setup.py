from setuptools import setup

setup(
    name='infraredPumping',
    version='1.0',
    description='Calculate pumping of rotational states by rovibrational cascade using billions of spectroscopic transitions from HITRAN, LAMDA, and GEISA',
    packages = [
        'pumping', 'utils'
    ],
    install_requires = [
        'numpy', 'astropy', 'astroquery', 'pandas', 'radis'
    ]
)