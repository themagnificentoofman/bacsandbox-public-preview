import os
from pathlib import Path
import shlex
import sys
import tempfile

import numpy
import pybind11
from setuptools import Extension, find_packages, setup
from setuptools._distutils.errors import (
    CompileError,
    DistutilsExecError,
    DistutilsOptionError,
    DistutilsPlatformError,
    LinkError,
)
from setuptools.command.build_ext import build_ext


ROOT = Path(__file__).resolve().parent
OPENMP_ENVIRONMENT_VARIABLE = 'BACSANDBOX_OPENMP'
OPENMP_POLICIES = frozenset(('auto', 'on', 'off'))
OPENMP_PROBE_SOURCE = r'''
#include <omp.h>

int main() {
    int thread_count = 0;
#pragma omp parallel reduction(+:thread_count)
    thread_count += 1;
    return thread_count > 0 && omp_get_max_threads() > 0 ? 0 : 1;
}
'''


class KernelBuildExt(build_ext):
    def build_extensions(self):
        policy = os.environ.get(OPENMP_ENVIRONMENT_VARIABLE, 'auto').strip().lower()
        if policy not in OPENMP_POLICIES:
            raise DistutilsOptionError(
                f'{OPENMP_ENVIRONMENT_VARIABLE} must be one of auto, on, or off; got {policy!r}'
            )

        if self.compiler.compiler_type == 'msvc':
            release_compile_args = ['/O2', '/fp:precise', '/std:c++17']
        else:
            release_compile_args = ['-O3', '-std=c++17']

        openmp_compile_args, openmp_link_args = self._select_openmp_flags(policy)
        openmp_enabled = bool(openmp_compile_args or openmp_link_args)
        for extension in self.extensions:
            extension.extra_compile_args = [
                *(getattr(extension, 'extra_compile_args', None) or []),
                *release_compile_args,
                *openmp_compile_args,
            ]
            extension.extra_link_args = [
                *(getattr(extension, 'extra_link_args', None) or []),
                *openmp_link_args,
            ]
            extension.define_macros = [
                *(getattr(extension, 'define_macros', None) or []),
                ('BACSANDBOX_OPENMP_ENABLED', '1' if openmp_enabled else '0'),
            ]

        super().build_extensions()

    def _select_openmp_flags(self, policy):
        if policy == 'off':
            self.announce('OpenMP probing disabled by BACSANDBOX_OPENMP=off.', level=2)
            return [], []

        failures = []
        for compile_args, link_args in self._openmp_candidates():
            try:
                self._probe_openmp(compile_args, link_args)
            except (CompileError, LinkError, DistutilsExecError) as exc:
                failures.append(f'compile={compile_args!r}, link={link_args!r}: {exc}')
                continue

            self.announce(
                f'OpenMP probe succeeded with compile={compile_args!r}, link={link_args!r}.',
                level=2,
            )
            return compile_args, link_args

        if policy == 'on':
            details = '\n'.join(failures) if failures else 'No OpenMP probe candidates were available.'
            raise DistutilsPlatformError(
                'BACSANDBOX_OPENMP=on requires OpenMP, but every compile/link probe failed:\n'
                + details
            )

        self.announce(
            'OpenMP was unavailable; building bacteria_c_kernel in serial mode.',
            level=3,
        )
        return [], []

    def _openmp_candidates(self):
        if self.compiler.compiler_type == 'msvc':
            return [(['/openmp'], [])]

        if sys.platform != 'darwin':
            return [(['-fopenmp'], ['-fopenmp'])]

        cppflags = shlex.split(os.environ.get('CPPFLAGS', ''))
        ldflags = shlex.split(os.environ.get('LDFLAGS', ''))
        candidates = []
        if cppflags or ldflags:
            candidates.extend((
                (
                    ['-Xpreprocessor', '-fopenmp', *cppflags],
                    [*ldflags, '-lomp'],
                ),
                (
                    ['-fopenmp', *cppflags],
                    ['-fopenmp', *ldflags],
                ),
            ))
        candidates.extend((
            (['-Xpreprocessor', '-fopenmp'], ['-lomp']),
            (
                ['-Xpreprocessor', '-fopenmp', '-I/opt/homebrew/opt/libomp/include'],
                ['-L/opt/homebrew/opt/libomp/lib', '-lomp'],
            ),
            (
                ['-Xpreprocessor', '-fopenmp', '-I/usr/local/opt/libomp/include'],
                ['-L/usr/local/opt/libomp/lib', '-lomp'],
            ),
        ))
        return candidates

    def _probe_openmp(self, compile_args, link_args):
        with tempfile.TemporaryDirectory(prefix='bacsandbox-openmp-') as temporary_directory:
            probe_directory = Path(temporary_directory)
            source_path = probe_directory / 'openmp_probe.cpp'
            source_path.write_text(OPENMP_PROBE_SOURCE, encoding='ascii')
            objects = self.compiler.compile(
                [str(source_path)],
                output_dir=str(probe_directory),
                extra_postargs=compile_args,
            )
            self.compiler.link_executable(
                objects,
                str(probe_directory / 'openmp_probe'),
                extra_postargs=link_args,
                target_lang='c++',
            )


setup(
    name='bacsandbox',
    version='0.1.0',
    description='Interactive bacterial ecosystem sandbox',
    long_description=(ROOT / 'README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    python_requires='>=3.9,<3.14',
    install_requires=[
        'pygame==2.6.1',
        'numpy>=1.24',
        'pandas>=2.0',
        'matplotlib>=3.7',
    ],
    packages=find_packages(exclude=('tests', 'tests.*')),
    include_package_data=False,
    package_data={'bacteria': ['models/*.json']},
    ext_modules=[
        Extension(
            'bacteria_c_kernel',
            sources=[
                'bacteria/src/ecm_kernel.cpp',
                'bacteria/src/metabolic_kernel.cpp',
                'bacteria/src/physics_kernel.cpp',
            ],
            include_dirs=[pybind11.get_include(), numpy.get_include()],
            depends=[
                'bacteria/src/ecm_kernel.h',
                'bacteria/src/metabolic_kernel.h',
                'bacteria/src/physics_kernel.h',
            ],
            language='c++',
        ),
    ],
    cmdclass={'build_ext': KernelBuildExt},
    zip_safe=False,
)