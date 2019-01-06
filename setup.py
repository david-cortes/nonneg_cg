try:
	from setuptools import setup
except:
	from distutils.core import setup
import numpy as np
from distutils.extension import Extension
from findblas.distutils import build_ext_with_blas

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
		build_ext_with_blas.build_extensions(self)

setup(
	name  = "nonnegcg",
	packages = ["nonnegcg"],
	version = '0.1',
	description = 'Conjugate-gradient optimizer subject to non-negativity constraints',
	author = 'David Cortes',
	author_email = 'david.cortes.rivera@gmail.com',
	url = 'https://github.com/david-cortes/nonneg_cg',
	keywords = ['optimization', 'conjugate gradient', 'polak-rubiere-polyak'],
	install_requires=[
		'numpy',
		'scipy',
		'cython',
		'findblas'
	],
	cmdclass = {'build_ext': build_ext_with_blas},
	data_files = [('include', ['include/nonnegcg.h', 'nonnegcg/nonnegcg.pxd'])],
	include_package_data = True,
	ext_modules = [
		Extension("nonnegcg._minimize", sources=["nonnegcg/pywrapper.pyx"], include_dirs=[np.get_include()], define_macros=[("_FOR_PYTON", None)]
			)]
	)
