from setuptools import Extension, setup

setup(ext_modules=[Extension("varquanms._native", ["src/varquanms/_native.c"],
                            libraries=["z"],
                            extra_compile_args=["-O3", "-fno-fast-math", "-ffp-contract=off"],
                            optional=True)])
