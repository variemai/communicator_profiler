# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install mpisee
#
# You can edit this file again by typing:
#
#     spack edit mpisee
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *
import os

class Mpisee(CMakePackage):
    """Mpisee: A communicator-centric profiler."""  # Description

    homepage = "https://github.com/variemai/communicator_profiler"
    git      = "https://github.com/variemai/communicator_profiler.git"
    maintainers = ['variemai']  # Add your GitHub username(s)

    version('3.1', branch='main')
    depends_on('mpi', type=('link'))
    depends_on('sqlite', type=('link'))

    # def cmake_args(self):
    #     args = [self.define_from_variant('NUM_BUCKETS', 'buckets')]
    #     args = [self.define_from_variant('BUCKETS', 'buckets')]
    #     return args

    # variant(
    #     "buckets",
    #     default="default",
    #     values=("default"),
    #     description="Number of buckets and the bucket sizes"
    # )

    # This install method is the same as the CMakeLists.txt
    def install(self, spec, prefix):
        mkdir(prefix.lib)
        mkdir(prefix.include)
        mkdir(prefix.bin)

        build_dir = self.build_directory
        seethrough_dir = join_path(self.stage.source_path, 'mpisee-through')

        # Install the mpisee-through.py script
        install(join_path(seethrough_dir, 'mpisee-through.py'), prefix.bin)
        # Make mpisee-through.py executable
        seethrough_script = join_path(prefix.bin, 'mpisee-through.py')
        os.chmod(seethrough_script, 0o755)  # Set permissions using os.chmod

        # Install the libmpisee.so library
        install(join_path(build_dir, 'libmpisee.so'), prefix.lib)
