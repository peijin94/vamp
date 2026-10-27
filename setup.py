from setuptools import Extension, setup

setup(ext_modules=[Extension("vamp._native", ["src/vamp/_native.c"],
                            libraries=["z"],
                            extra_compile_args=["-O3", "-fno-fast-math", "-ffp-contract=off"],
                            optional=True)])
