try:
	from setuptools import setup
	from setuptools import Extension
except:
	from distutils.core import setup
	from distutils.extension import Extension
import numpy as np, warnings
from findblas.distutils import build_ext_with_blas
from sys import platform

## https://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used
class build_ext_subclass( build_ext_with_blas ):
	def build_extensions(self):
		compiler = self.compiler.compiler_type
		if compiler == 'msvc': # visual studio
			for e in self.extensions:
				e.extra_compile_args += ['/O2', '/openmp']
		else: # everything else that cares about following standards
			for e in self.extensions:
				e.extra_compile_args += ['-O2', '-fopenmp', '-march=native', '-std=c99']
				e.extra_link_args += ['-fopenmp']

		if platform[:3] == "dar":
			apple_msg  = "\n\n\nMacOS detected. Package will be built without multi-threading capabilities, "
			apple_msg += "due to Apple's lack of OpenMP support in default Xcode installs. In order to enable it, "
			apple_msg += "install the package directly from GitHub: https://www.github.com/david-cortes/nonneg_cg\n"
			apple_msg += "And modify the setup.py file where this message is shown. "
			apple_msg += "You'll also need an OpenMP-capable compiler.\n\n\n"
			warnings.warn(apple_msg)
			for e in self.extensions:
				e.extra_compile_args = [arg for arg in extra_compile_args if arg != '-fopenmp']
				e.extra_link_args	= [arg for arg in extra_link_args	if arg != '-fopenmp']


		build_ext_with_blas.build_extensions(self)

setup(
	name  = "nonnegcg",
	packages = ["nonnegcg"],
	version = '0.1.6',
	description = 'Conjugate-gradient optimizer subject to non-negativity constraints',
	author = 'David Cortes',
	author_email = 'david.cortes.rivera@gmail.com',
	url = 'https://github.com/david-cortes/nonneg_cg',
	keywords = ['optimization', 'conjugate gradient', 'polak-ribiere-polyak'],
	install_requires=[
		'numpy',
		'scipy',
		'cython',
		'findblas'
	],
	cmdclass = {'build_ext': build_ext_subclass},
	data_files = [('include', ['include/nonnegcg.h', 'nonnegcg/nonnegcg.pxd'])],
	include_package_data = True,
	ext_modules = [
		Extension("nonnegcg._minimize", sources=["nonnegcg/pywrapper.pyx"], include_dirs=[np.get_include()], define_macros=[("_FOR_PYTHON", None)]
			)]
	)
